/*
 * Copyright 2022 - UPMEM
 */

#include <cassert>
#include <cstdio>
#include <fstream>

#include "dpu_common.hpp"

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

std::vector<NW_dpu_input> fair_dispatch(const Sets &data, size_t nr_of_dpus)
{
    printf("\n> Using fair dispatch:\n");
    printf("   sets number: %zu\n", data.size());
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

    printf("   margin: %.1f%%\n", static_cast<float>(margin * 100) / static_cast<float>(mean_compute_load));

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

    printf("   dpu not used: %zu\n", nr_of_dpus - 1 - dpu_index);

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