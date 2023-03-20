/*
 * Copyright 2022 - UPMEM
 */

#ifndef E6039E80_5D9F_462C_ACAE_D977B65797AC
#define E6039E80_5D9F_462C_ACAE_D977B65797AC

#include "../../src/types.hpp"

/**
 * @brief Structure regrouping all data to be send to a dpu
 *
 */
typedef struct NW_dpu_input
{
    NwMetadataDPU metadata{};      /// metadata describing the sequences
    CompressedSequences sequences{};       /// all sequences in one single buffer
    std::vector<uint32_t> cigar_indexes{}; /// start index of each cigar memory space
} NW_dpu_input;

typedef struct NW_score_input
{
    NwMetadataDPU metadata{};
    CompressedSequences sequences{};
} NW_score_input;

std::vector<nw_t> dpu_pipeline(std::string dpu_bin_path, const NW_Parameters &p, size_t nr_dpu, const Sets &sets);
std::vector<int> dpu_pipeline_16s(std::string dpu_bin_path, const NW_Parameters &params, size_t ndpu, const Set &set);

#endif /* E6039E80_5D9F_462C_ACAE_D977B65797AC */
