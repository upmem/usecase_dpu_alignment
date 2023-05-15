/*
 * Copyright 2022 - UPMEM
 */

#include <algorithm>
#include <cassert>
#include <cstdio>
#include <fstream>

#include "dpu_common.hpp"

extern "C"
{
#include <dpu.h>
}

uint32_t compressed_size(const Sequence &seq)
{
    static constexpr size_t alignment = 8;
    const auto compressed_size = (seq.size() + 3) / 4;
    if (compressed_size % 8 == 0)
        return compressed_size;
    return compressed_size + (alignment - (compressed_size % alignment));
}

CompressedSequence compress_sequence(const Sequence &seq)
{
    uint32_t csize = compressed_size(seq);
    CompressedSequence cseq(csize);

    for (size_t i = 0; i < cseq.size(); i++)
    {
        size_t seq_id = i * 4;
        uint8_t c4n = seq[seq_id];
        c4n |= (seq[seq_id + 1] << 2);
        c4n |= (seq[seq_id + 2] << 4);
        c4n |= (seq[seq_id + 3] << 6);

        cseq[i] = c4n;
    }

    return cseq;
}

CompressedSet compress_set(const Set &set)
{
    CompressedSet cset(set.size());

#pragma omp parallel for
    for (size_t i = 0; i < set.size(); i++)
        cset[i] = compress_sequence(set[i]);

    return cset;
}

void push_back(CompressedSequences &cseqs, const CompressedSequence &cseq)
{
    size_t begin = cseqs.size();
    cseqs.resize(begin + cseq.size());

    for (size_t i = 0; i < cseq.size(); i++)
        cseqs[begin + i] = cseq[i];
}

inline auto sorted_map(const Sets &data)
{
    std::vector<std::tuple<int, size_t, size_t>> index{data.size()};

    for (size_t i = 0; i < index.size(); i++)
        index[i] = {i, count_compute_load(data[i]), 0};

    std::sort(index.begin(), index.end(), [](const auto &a, const auto &b)
              { return std::get<1>(a) > std::get<1>(b); });

    return index;
}

auto bucket_sets(const Sets &data, auto &index, size_t n)
{
    std::vector<Sets> dpu_sets(n);
    std::vector<size_t> dpu_loads(n);

    for (auto &[i, load, d] : index)
    {
        auto min_index = std::distance(dpu_loads.begin(), std::min_element(dpu_loads.begin(), dpu_loads.end()));
        dpu_sets[min_index].push_back(data[i]);
        dpu_loads[min_index] += load;
        d = min_index;
    }

    return dpu_sets;
}

inline void cpu_to_dpu(const Sets &sets, NwInputCigar &dpu_input)
{
    assert(sets.size() <= SCORE_METADATA_MAX_NUMBER_OF_SET &&
           "Too many sets for DPU!\n");

    auto &meta = dpu_input.metadata;
    meta.number_of_sets = static_cast<uint32_t>(sets.size());

    uint32_t idx = 0;
    uint32_t seq_idx = 0;
    uint32_t cigar_index = 0;
    uint32_t cigar_offset = 0;

    for (size_t set_id = 0; set_id < sets.size(); set_id++)
        meta.set_sizes[set_id] = static_cast<uint8_t>(sets[set_id].size());

    for (const auto &set : sets)
    {
        for (const auto &seq : set)
        {
            auto cseq = compress_sequence(seq);
            meta.lengths[seq_idx] = static_cast<uint16_t>(seq.size());
            meta.indexes[seq_idx++] = idx;
            idx += cseq.size();
            push_back(dpu_input.sequences, cseq);
        }
        assert(dpu_input.sequences.size() < SCORE_MAX_SEQUENCES_TOTAL_SIZE &&
               "dpu sequence buffer overflow!\n");

        for (size_t i = 0; i < set.size() - 1; i++)
            for (size_t j = i + 1; j < set.size(); j++)
            {
                dpu_input.cigar_indexes[cigar_offset++] = cigar_index;
                auto max_cigar_size = set[i].size() + set[j].size();

                assert(set[i].size() + set[j].size() < UINT16_MAX && "cigar is to big for uint16_t\n");

                max_cigar_size += (8 - (max_cigar_size % 8));
                cigar_index += max_cigar_size;

                assert(cigar_index < MAX_CIGAR_SIZE && "not enough space for cigar on dpu!\n");
                assert(cigar_offset <= METADATA_MAX_NUMBER_OF_SCORES && "too much cigar to compute!\n");
            }
    }
}

auto dispatchv2(const Sets &data, size_t nr_of_dpus, const NwParameters &p, std::vector<std::tuple<int, size_t, size_t>> &index)
{
    index = sorted_map(data);
    auto dpu_sets = bucket_sets(data, index, nr_of_dpus);

    std::vector<NwInputCigar> dpu_input(nr_of_dpus);
#pragma omp parallel for num_threads(2)
    for (size_t i = 0; i < nr_of_dpus; i++)
    {
        dpu_input[i].metadata.match = p.match;
        dpu_input[i].metadata.mismatch = p.mismatch;
        dpu_input[i].metadata.gap_opening = p.gap_opening;
        dpu_input[i].metadata.gap_extension = p.gap_extension;
        dpu_input[i].sequences.reserve(SCORE_MAX_SEQUENCES_TOTAL_SIZE);
        dpu_input[i].cigar_indexes.resize(METADATA_MAX_NUMBER_OF_SCORES);

        cpu_to_dpu(dpu_sets[i], dpu_input[i]);
    }

    return dpu_input;
}

auto dpu_to_cpu(std::vector<std::vector<NwType>> &res, const NwInputCigar &input, const NwCigarOutput &output, size_t mi)
{
    const auto &meta = input.metadata;
    res.resize(meta.number_of_sets);

    auto comb = [](auto n)
    { return n * (n - 1) / 2; };

    int dpu_offset = 0;

    for (size_t set_id = 0; set_id < meta.number_of_sets; set_id++)
    {
        int set_offset = 0;
        res[set_id].resize(comb(meta.set_sizes[set_id]));
        for (int i = 0; i < meta.set_sizes[set_id]; i++)
            for (int j = i + 1; j < meta.set_sizes[set_id]; j++)
            {
                res[set_id][set_offset].score = output.scores[dpu_offset];
                res[set_id][set_offset].dpu_offset = dpu_offset;
                res[set_id][set_offset].mi = mi;
                dpu_offset++;
                set_offset++;
            }
    }

    return res;
}

std::vector<NwType> dpus_to_cpu(const auto &sets, const auto &inputs, const auto &outputs, const auto &cigars, const std::vector<std::tuple<int, size_t, size_t>> &index)
{
    std::vector<std::vector<std::vector<NwType>>> dpu_res(inputs.size());
#pragma omp parallel for num_threads(2)
    for (size_t i = 0; i < inputs.size(); i++)
        dpu_to_cpu(dpu_res[i], inputs[i], outputs[i], i);

    std::vector<std::vector<NwType>> cpu_res{sets.size()};
    for (const auto &[i, load, d] : index)
    {
        cpu_res[i] = dpu_res[d].front();
        dpu_res[d].erase(dpu_res[d].begin());
    }

    std::vector<NwType> cpu_output(count_unique_pair(sets));
    size_t g = 0;
    for (const auto &res : cpu_res)
        for (const auto &r : res)
        {
            cpu_output[g] = r;
            cpu_output[g].cigar.resize(outputs[r.mi].lengths[r.dpu_offset]);
            cpu_output[g].cigar.assign(&cigars[r.mi][inputs[r.mi].cigar_indexes[r.dpu_offset]], outputs[r.mi].lengths[r.dpu_offset]);
            g++;
        }

    return cpu_output;
}

void send_cigar_input(dpu_set_t &dpu_set, std::vector<NwInputCigar> &inputs)
{
    dpu_set_t dpu{};
    uint32_t each_dpu = 0;

    DPU_FOREACH(dpu_set, dpu, each_dpu)
    {
        DPU_ASSERT(dpu_prepare_xfer(dpu, inputs[each_dpu].sequences.data()));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_TO_DPU, "sequences", 0,
                             SCORE_MAX_SEQUENCES_TOTAL_SIZE, DPU_XFER_ASYNC));

    DPU_FOREACH(dpu_set, dpu, each_dpu)
    {
        DPU_ASSERT(dpu_prepare_xfer(dpu, &inputs[each_dpu].metadata));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_TO_DPU, "metadata", 0,
                             sizeof(NwMetadataDPU), DPU_XFER_ASYNC));

    DPU_FOREACH(dpu_set, dpu, each_dpu)
    {
        DPU_ASSERT(dpu_prepare_xfer(dpu, inputs[each_dpu].cigar_indexes.data()));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_TO_DPU, "cigar_indexes", 0,
                             METADATA_MAX_NUMBER_OF_SCORES * sizeof(uint32_t), DPU_XFER_ASYNC));
}

void gather_cigar_output(dpu_set_t &dpu_set, std::vector<NwCigarOutput> &outputs, std::vector<std::vector<char>> &cigars)
{
    dpu_set_t dpu{};
    uint32_t each_dpu = 0;

    for (auto &cigar : cigars)
        cigar.resize(MAX_CIGAR_SIZE);

    DPU_FOREACH(dpu_set, dpu, each_dpu)
    {
        DPU_ASSERT(dpu_prepare_xfer(dpu, &outputs[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_FROM_DPU, "output", 0,
                             sizeof(NwCigarOutput), DPU_XFER_ASYNC));
    DPU_FOREACH(dpu_set, dpu, each_dpu)
    {
        DPU_ASSERT(dpu_prepare_xfer(dpu, cigars[each_dpu].data()));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_FROM_DPU, "cigars", 0,
                             MAX_CIGAR_SIZE, DPU_XFER_ASYNC));
}

dpu_set_t init_dpu(const std::filesystem::path &filename, size_t count)
{
    dpu_set_t dpus{};
    if (not std::filesystem::exists(filename))
        exit("File " + filename.string() + " not found.");

    DPU_ASSERT(dpu_alloc(count, NULL, &dpus));
    DPU_ASSERT(dpu_load(dpus, filename.c_str(), NULL));
    return dpus;
}

std::vector<NwType> dpu_cigar_pipeline(std::string dpu_bin_path, const NwParameters &p, size_t nr_dpu, const Sets &sets)
{
    Timer dispatch_time{};
    std::vector<std::tuple<int, size_t, size_t>> index;
    auto dpu_inputs = dispatchv2(sets, nr_dpu, p, index);
    dispatch_time.Print("  ");

    std::vector<NwCigarOutput> dpu_outputs(nr_dpu);
    std::vector<std::vector<char>> dpu_cigars(nr_dpu);

    printf("\nAllocating DPUs\n");
    auto dpus = init_dpu(dpu_bin_path, nr_dpu);

    printf("Send, Launch, Gather (Async)\n");
    send_cigar_input(dpus, dpu_inputs);
    DPU_ASSERT(dpu_launch(dpus, DPU_ASYNCHRONOUS));
    gather_cigar_output(dpus, dpu_outputs, dpu_cigars);
    dpu_sync(dpus);

    // Add to dump DPU counters and analyse their individual workload.
    // dump_to_file("counters.txt", dpu_outputs, [](const auto &e) { return e.perf_counter; });

    return dpus_to_cpu(sets, dpu_inputs, dpu_outputs, dpu_cigars, index);
}

auto Set_to_dpuSet(const Set &data, const NwParameters &params)
{
    NwInputScore dpu_input;

    dpu_input.sequences.reserve(SCORE_MAX_SEQUENCES_TOTAL_SIZE);
    dpu_input.metadata.match = params.match;
    dpu_input.metadata.mismatch = params.mismatch;
    dpu_input.metadata.gap_extension = params.gap_extension;
    dpu_input.metadata.gap_opening = params.gap_opening;

    size_t dpu_index = 0;
    size_t seq_id = 0;

    auto cset = compress_set(data);

    for (size_t i = 0; i < cset.size(); i++)
    {
        assert(seq_id < DPU_MAX_NUMBER_OF_SEQUENCES_MRAM && "Set is too big!\n");
        dpu_input.sequence_metadata.lengths[seq_id] = data[i].size();
        dpu_input.sequence_metadata.indexes[seq_id] = dpu_index;
        dpu_index += cset[i].size();
        seq_id++;
    }

    for (const auto &seq : cset)
        push_back(dpu_input.sequences, seq);

    assert(dpu_input.sequences.size() < SCORE_MAX_SEQUENCES_TOTAL_SIZE &&
           "dpu sequence buffer overflow!\n");

    return dpu_input;
}

void send_input_16s(dpu_set_t &dpu_set, NwInputScore &inputs, std::vector<ComparisonMetadata> &meta)
{

    DPU_ASSERT(dpu_broadcast_to(dpu_set, "sequences", 0, inputs.sequences.data(), inputs.sequences.size(), DPU_XFER_ASYNC));
    DPU_ASSERT(dpu_broadcast_to(dpu_set, "metadata", 0, &(inputs.metadata), sizeof(NwMetadataDPU), DPU_XFER_ASYNC));
    DPU_ASSERT(dpu_broadcast_to(dpu_set, "sequence_metadata", 0, &(inputs.sequence_metadata), sizeof(NwSequenceMetadataMram), DPU_XFER_ASYNC));

    dpu_set_t dpu{};
    uint32_t each_dpu = 0;

    DPU_FOREACH(dpu_set, dpu, each_dpu)
    {
        DPU_ASSERT(dpu_prepare_xfer(dpu, &(meta[each_dpu])));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_TO_DPU, "meta_index", 0,
                             sizeof(ComparisonMetadata), DPU_XFER_ASYNC));
}

void gather_output_16s(dpu_set_t &dpu_set, std::vector<NwScoreOutput> &outputs)
{
    dpu_set_t dpu{};
    uint32_t each_dpu = 0;

    DPU_FOREACH(dpu_set, dpu, each_dpu)
    {
        DPU_ASSERT(dpu_prepare_xfer(dpu, &outputs[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_FROM_DPU, "output", 0,
                             sizeof(NwScoreOutput), DPU_XFER_ASYNC));
}

void update_meta(ComparisonMetadata &meta, int &rest)
{
    for (uint32_t i = 0; i < meta.count; i++)
    {
        meta.start_col++;
        if (meta.start_col >= meta.size)
        {
            meta.start_row++;
            meta.start_col = meta.start_row + 1;
        }
    }

    if (rest-- == 0)
        meta.count--;
}

auto dispatch(const Set &set, size_t nr_dpu)
{
    auto total_size = sum_integers(set.size());
    auto mean = total_size / nr_dpu;
    int rest = static_cast<int>(total_size % nr_dpu);

    std::vector<ComparisonMetadata> dpu_metadata(nr_dpu);

    printf("  mean: %lu\n", mean);

    ComparisonMetadata meta{
        0,
        1,
        static_cast<uint32_t>(mean + (rest-- != 0 ? 1 : 0)),
        static_cast<uint32_t>(set.size())};

    for (auto &dpu : dpu_metadata)
    {
        dpu = meta;
        update_meta(meta, rest);
    }

    return dpu_metadata;
}

std::vector<int> dpu_16s_pipeline(std::string dpu_bin_path, const NwParameters &params, size_t ndpu, const Set &set)
{
    auto metadata = dispatch(set, ndpu);
    auto compressed_set = Set_to_dpuSet(set, params);

    printf("\nAllocating DPUs\n");
    auto dpus = init_dpu(dpu_bin_path, ndpu);

    printf("Send, Launch, Gather (Async)\n");
    send_input_16s(dpus, compressed_set, metadata);
    DPU_ASSERT(dpu_launch(dpus, DPU_ASYNCHRONOUS));

    std::vector<NwScoreOutput> outputs(ndpu);
    gather_output_16s(dpus, outputs);

    std::vector<int> results_dpu{};
    results_dpu.reserve(sum_integers(set.size()));
    dpu_sync(dpus);

    for (const auto &dpu : outputs)
        for (size_t i = 0; i < dpu.nr_score; i++)
            results_dpu.push_back(dpu.scores[i]);

    // Add to dump DPU counters and analyse their individual workload.
    dump_to_file("counters.txt", outputs, [](const auto &e)
                 { return e.perf_counter; });

    return results_dpu;
}
