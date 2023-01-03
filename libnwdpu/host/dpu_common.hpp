/*
 * Copyright 2022 - UPMEM
 */

#ifndef E6039E80_5D9F_462C_ACAE_D977B65797AC
#define E6039E80_5D9F_462C_ACAE_D977B65797AC

extern "C"
{
#include <dpu.h>
}

#include <algorithm>
#include <fstream>

#include "../../src/types.hpp"

/**
 * @brief Structure regrouping all data to be send to a dpu
 *
 */
typedef struct NW_dpu_input
{
    NW_dpu_metadata_input metadata{};      /// metadata describing the sequences
    CompressedSequences sequences{};       /// all sequences in one single buffer
    std::vector<uint32_t> cigar_indexes{}; /// start index of each cigar memory space
} NW_dpu_input;

/**
 * @brief Returns size of compressed sequence with 8-aligned padding
 *
 * @param seq
 * @return uint32_t
 */
uint32_t compressed_size(const Sequence &seq);

uint32_t add_sequence(std::vector<uint8_t> &dseq, const Sequence &seq, uint32_t idx);

std::vector<NW_dpu_input> fair_dispatch(const Sets &data, size_t nr_of_dpus);

void send_cigar_input(dpu_set_t &dpu_set, std::vector<NW_dpu_input> &inputs);
void gather_cigar_output(dpu_set_t &dpu_set, std::vector<NW_dpu_output> &outputs, std::vector<std::vector<uint8_t>> &cigars);

static inline auto dpu_pipeline(std::string dpu_bin_path, const NW_Parameters &p, size_t nr_dpu, const Sets &sets)
{

    dpu_set_t dpus{};
    DPU_ASSERT(dpu_alloc(nr_dpu, NULL, &dpus));
    DPU_ASSERT(dpu_load(dpus, dpu_bin_path.c_str(), NULL));

    Timer dispatch_time{};
    auto dpu_inputs = fair_dispatch(sets, nr_dpu);
    printf("dispatch:\n");
    dispatch_time.print();

    for (auto &dpu : dpu_inputs)
    {
        dpu.metadata.match = p.match;
        dpu.metadata.mismatch = p.mismatch;
        dpu.metadata.gap_opening = p.gap_opening1;
        dpu.metadata.gap_extension = p.gap_extension1;
    }

    printf("> transfert\n");
    send_cigar_input(dpus, dpu_inputs);

    printf("> Launch dpu\n");
    DPU_ASSERT(dpu_launch(dpus, DPU_ASYNCHRONOUS));

    std::vector<NW_dpu_output> outputs(nr_dpu);
    std::vector<std::vector<uint8_t>> dpu_cigars(nr_dpu);
    gather_cigar_output(dpus, outputs, dpu_cigars);

    dpu_sync(dpus);

    printf("   done!\n");
    fflush(stdout);

    printf("perf counter: %lu\n", outputs[0].perf_counter);

    std::ofstream fcounter("counter.txt");

    for (const auto &e : outputs)
        fcounter << e.perf_counter << '\n';

    std::vector<nw_t> res(count_unique_pair(sets));

    auto comb = [](auto n)
    { return n * (n - 1) / 2; };

    size_t i = 0;
    for (size_t d = 0; d < dpu_inputs.size(); d++)
    {
        size_t s = 0;
        for (size_t n = 0; n < dpu_inputs[d].metadata.number_of_sets; n++)
            for (int ns = 0; ns < comb(dpu_inputs[d].metadata.set_sizes[n]); ns++)
            {
                res[i].score = outputs[d].scores[s];
                res[i].cigar.resize(outputs[d].lengths[s]);
                res[i].cigar.assign(reinterpret_cast<char *>(&dpu_cigars[d][dpu_inputs[d].cigar_indexes[s]]), outputs[d].lengths[s]);

                s++;
                i++;
            }
    }

    return res;
}

#endif /* E6039E80_5D9F_462C_ACAE_D977B65797AC */
