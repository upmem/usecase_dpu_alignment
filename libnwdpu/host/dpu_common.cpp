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

inline bool load_can_fit(const Sets &sets, size_t mean_compute_load, size_t nr_of_dpus, size_t margin)
{
    size_t dpu_index = 0;
    std::vector<size_t> dpu_compute_load(nr_of_dpus, 0);
    std::vector<size_t> set_compute_load{};

    for (const auto &set : sets)
        set_compute_load.emplace_back(count_compute_load(set));

    for (size_t i = 0; i < sets.size(); i++)
    {
        if (dpu_compute_load[dpu_index] + set_compute_load[i] > (mean_compute_load + margin))
            dpu_index++;

        if (dpu_index >= nr_of_dpus)
            return false;

        dpu_compute_load[dpu_index] += set_compute_load[i];
    }

    return true;
}

inline size_t find_lowest_limit(const Sets &sets, size_t mean_compute_load, size_t nr_of_dpus)
{
    uint32_t mil = mean_compute_load / 1000;
    for (size_t margin = 0;; margin += mil)
        if (load_can_fit(sets, mean_compute_load, nr_of_dpus, margin))
            return margin;
}

/**
 * @brief check if current dpu can fit the set compute load.
 * Go to the next dpu if not. Assert if not enough DPUs.
 *
 * @param set set to add
 * @param dpu_compute_load current load of DPUs.
 * @param mean mean load of all set per DPU to target
 * @param margin smallest working margin
 * @param idx current DPU index
 * @param nr_dpu total number of DPUs
 * @return size_t
 */
static inline size_t add_dpu_load(const auto &set, auto &dpu_compute_load, size_t mean, size_t margin, size_t idx, size_t nr_dpu)
{
    auto set_compute_load = count_compute_load(set);

    if (dpu_compute_load[idx] + set_compute_load > (mean + margin))
        idx++;

    assert(idx < nr_dpu && "Could not fit the dataset in dpus fairly, not enough dpus!\n");

    dpu_compute_load[idx] += set_compute_load;

    return idx;
}

std::vector<NW_dpu_input> fair_dispatch(const Sets &data, size_t nr_of_dpus, const NW_Parameters &p)
{
    printf("\nDispatch:\n");
    std::vector<NW_dpu_input> dpu_input(nr_of_dpus);
    for (auto &d : dpu_input)
    {
        d.sequences.reserve(SCORE_MAX_SEQUENCES_TOTAL_SIZE);
        d.cigar_indexes.resize(METADATA_MAX_NUMBER_OF_SCORES);
    }

    assert(data.size() <= nr_of_dpus * SCORE_METADATA_MAX_NUMBER_OF_SET &&
           "Too many sets per DPU!\n");

    auto mean_compute_load = count_compute_load(data) / nr_of_dpus;
    size_t dpu_index = 0;
    std::vector<size_t> dpu_offset(nr_of_dpus, 0);
    std::vector<size_t> dpu_indexes(nr_of_dpus, 0);
    std::vector<size_t> dpu_cigar_indexes_offset(nr_of_dpus, 0);
    std::vector<size_t> dpu_current_cigar_index(nr_of_dpus, 0);
    std::vector<size_t> dpu_compute_load(nr_of_dpus, 0);

    size_t margin = find_lowest_limit(data, mean_compute_load, nr_of_dpus);

    printf("  margin: %.1f%%\n", static_cast<float>(margin * 100) / static_cast<float>(mean_compute_load));

    for (const auto &set : data)
    {
        auto cset = compress_set(set);

        dpu_index = add_dpu_load(set, dpu_compute_load, mean_compute_load, margin, dpu_index, nr_of_dpus);

        auto &meta = dpu_input[dpu_index].metadata;

        assert(meta.number_of_sets < SCORE_METADATA_MAX_NUMBER_OF_SET &&
               "Too much sets in one dpu!\n");

        auto off = dpu_offset[dpu_index];
        for (size_t i = 0; i < set.size(); i++)
        {
            assert(off + i < DPU_MAX_NUMBER_OF_SEQUENCES &&
                   "Too much sequence in one dpu!\n");
            meta.lengths[off + i] = set[i].size();
            meta.indexes[off + i] = dpu_indexes[dpu_index];
            dpu_indexes[dpu_index] += cset[i].size();
        }
        dpu_offset[dpu_index] += set.size();

        for (size_t i = 0; i < set.size() - 1; i++)
            for (size_t j = i + 1; j < set.size(); j++)
            {
                dpu_input[dpu_index].cigar_indexes[dpu_cigar_indexes_offset[dpu_index]++] = dpu_current_cigar_index[dpu_index];
                uint32_t max_cigar_size = set[i].size() + set[j].size();

                assert(set[i].size() + set[j].size() < UINT16_MAX && "cigar is to big for uint16_t\n");

                // max_cigar_size -= max_cigar_size / 4;
                max_cigar_size += (8 - (max_cigar_size % 8));
                dpu_current_cigar_index[dpu_index] += max_cigar_size;

                assert(dpu_current_cigar_index[dpu_index] < MAX_CIGAR_SIZE && "not enough space for cigar on dpu!\n");
                assert(dpu_cigar_indexes_offset[dpu_index] <= METADATA_MAX_NUMBER_OF_SCORES && "too much cigar to compute!\n");
            }

        meta.set_sizes[meta.number_of_sets] = set.size();
        meta.number_of_sets++;

        for (const auto &cseq : cset)
            push_back(dpu_input[dpu_index].sequences, cseq);

        assert(dpu_input[dpu_index].sequences.size() < SCORE_MAX_SEQUENCES_TOTAL_SIZE &&
               "dpu sequence buffer overflow!\n");
    }

    printf("  dpu not used: %zu\n", nr_of_dpus - 1 - dpu_index);

    for (auto &dpu : dpu_input)
    {
        dpu.metadata.match = p.match;
        dpu.metadata.mismatch = p.mismatch;
        dpu.metadata.gap_opening = p.gap_opening;
        dpu.metadata.gap_extension = p.gap_extension;
    }

    return dpu_input;
}

void send_cigar_input(dpu_set_t &dpu_set, std::vector<NW_dpu_input> &inputs)
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
                             sizeof(NW_dpu_metadata_input), DPU_XFER_ASYNC));

    DPU_FOREACH(dpu_set, dpu, each_dpu)
    {
        DPU_ASSERT(dpu_prepare_xfer(dpu, inputs[each_dpu].cigar_indexes.data()));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_TO_DPU, "cigar_indexes", 0,
                             METADATA_MAX_NUMBER_OF_SCORES * sizeof(uint32_t), DPU_XFER_ASYNC));
}

void gather_cigar_output(dpu_set_t &dpu_set, std::vector<NW_dpu_output> &outputs, std::vector<std::vector<char>> &cigars)
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
                             sizeof(NW_dpu_output), DPU_XFER_ASYNC));
    DPU_FOREACH(dpu_set, dpu, each_dpu)
    {
        DPU_ASSERT(dpu_prepare_xfer(dpu, cigars[each_dpu].data()));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_FROM_DPU, "cigars", 0,
                             MAX_CIGAR_SIZE, DPU_XFER_ASYNC));
}

std::vector<nw_t> dpu_to_cpu_format(const auto &sets, const auto &inputs, const auto &outputs, const auto &cigars)
{
    std::vector<nw_t> cpu_output(count_unique_pair(sets));

    auto comb = [](auto n)
    { return n * (n - 1) / 2; };

    size_t i = 0;
    for (size_t d = 0; d < inputs.size(); d++)
    {
        size_t s = 0;
        for (size_t n = 0; n < inputs[d].metadata.number_of_sets; n++)
            for (int ns = 0; ns < comb(inputs[d].metadata.set_sizes[n]); ns++)
            {
                cpu_output[i].score = outputs[d].scores[s];
                cpu_output[i].cigar.resize(outputs[d].lengths[s]);
                cpu_output[i].cigar.assign(&cigars[d][inputs[d].cigar_indexes[s]], outputs[d].lengths[s]);

                s++;
                i++;
            }
    }

    return cpu_output;
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

std::vector<nw_t> dpu_pipeline(std::string dpu_bin_path, const NW_Parameters &p, size_t nr_dpu, const Sets &sets)
{
    Timer dispatch_time{};
    auto dpu_inputs = fair_dispatch(sets, nr_dpu, p);
    dispatch_time.print("  ");

    std::vector<NW_dpu_output> dpu_outputs(nr_dpu);
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

    return dpu_to_cpu_format(sets, dpu_inputs, dpu_outputs, dpu_cigars);
}

auto Set_to_dpuSet(const Set &data, const NW_Parameters &params)
{
    NW_score_input dpu_input;

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
        assert(seq_id < 1024 && "Set is too big!\n");
        dpu_input.metadata.lengths[seq_id] = data[i].size();
        dpu_input.metadata.indexes[seq_id] = dpu_index;
        dpu_index += cset[i].size();
        seq_id++;
    }

    for (const auto &seq : cset)
        push_back(dpu_input.sequences, seq);

    assert(dpu_input.sequences.size() < SCORE_MAX_SEQUENCES_TOTAL_SIZE &&
           "dpu sequence buffer overflow!\n");

    return dpu_input;
}

void send_input_16s(dpu_set_t &dpu_set, NW_score_input &inputs, std::vector<Metadata_index> &meta)
{

    DPU_ASSERT(dpu_broadcast_to(dpu_set, "sequences", 0, inputs.sequences.data(), inputs.sequences.size(), DPU_XFER_ASYNC));
    DPU_ASSERT(dpu_broadcast_to(dpu_set, "metadata", 0, &(inputs.metadata), sizeof(NW_dpu_metadata_input), DPU_XFER_ASYNC));

    dpu_set_t dpu{};
    uint32_t each_dpu = 0;

    DPU_FOREACH(dpu_set, dpu, each_dpu)
    {
        DPU_ASSERT(dpu_prepare_xfer(dpu, &(meta[each_dpu])));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_TO_DPU, "meta_index", 0,
                             sizeof(Metadata_index), DPU_XFER_ASYNC));
}

void gather_output_16s(dpu_set_t &dpu_set, std::vector<NW_score_output> &outputs)
{
    dpu_set_t dpu{};
    uint32_t each_dpu = 0;

    DPU_FOREACH(dpu_set, dpu, each_dpu)
    {
        DPU_ASSERT(dpu_prepare_xfer(dpu, &outputs[each_dpu]));
    }
    DPU_ASSERT(dpu_push_xfer(dpu_set, DPU_XFER_FROM_DPU, "output", 0,
                             sizeof(NW_score_output), DPU_XFER_ASYNC));
}

void update_meta(Metadata_index &meta, int &rest)
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

    std::vector<Metadata_index> dpu_metadata(nr_dpu);

    Metadata_index meta{
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

std::vector<int> dpu_pipeline_16s(std::string dpu_bin_path, const NW_Parameters &params, size_t ndpu, const Set &set)
{
    auto metadata = dispatch(set, ndpu);
    auto compressed_set = Set_to_dpuSet(set, params);

    printf("\nAllocating DPUs\n");
    auto dpus = init_dpu(dpu_bin_path, ndpu);

    printf("Send, Launch, Gather (Async)\n");
    send_input_16s(dpus, compressed_set, metadata);
    DPU_ASSERT(dpu_launch(dpus, DPU_ASYNCHRONOUS));

    std::vector<NW_score_output> outputs(ndpu);
    gather_output_16s(dpus, outputs);

    std::vector<int> results_dpu{};
    results_dpu.reserve(sum_integers(set.size()));
    dpu_sync(dpus);

    for (const auto &dpu : outputs)
        for (size_t i = 0; i < dpu.nr_score; i++)
            results_dpu.push_back(dpu.scores[i]);

    // Add to dump DPU counters and analyse their individual workload.
    // dump_to_file("counters.txt", outputs, [](const auto &e) { return e.perf_counter; });

    return results_dpu;
}