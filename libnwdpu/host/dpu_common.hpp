/*
 * Copyright 2022 - UPMEM
 */

#ifndef E6039E80_5D9F_462C_ACAE_D977B65797AC
#define E6039E80_5D9F_462C_ACAE_D977B65797AC

#include "../../src/types.hpp"

/**
 * @brief Structure regrouping all data to be send to a dpu for Cigar
 *
 */
struct NwInputCigar
{
    /// @brief Input data for Cigar NW
    NwMetadataDPU metadata{};              /// metadata describing the sequences
    CompressedSequences sequences{};       /// all sequences in are compressed into one single buffer
    std::vector<uint32_t> cigar_indexes{}; /// start index of each cigar memory space
};

/**
 * @brief Structure regrouping all data to be send to a dpu for Score
 *
 */
struct NwInputScore
{
    /// @brief Input data for Score NW
    NwMetadataDPU metadata{};        /// metadata describing the sequences
    CompressedSequences sequences{}; /// all sequences in are compressed into one single buffer
    NwSequenceMetadataMram sequence_metadata{};
};

/**
 * @brief DPU pipeline for CIGAR
 *
 * @param dpu_bin_path DPU binary path
 * @param params NW parameters
 * @param ranks Number of ranks to use
 * @param sets Dataset
 * @return std::vector<NwType>
 */
std::vector<NwType> dpu_cigar_pipeline(std::filesystem::path dpu_bin_path, const NwParameters &params, size_t ranks, const Sets &sets);

/**
 * @brief DPU pipeline for score
 *
 * @param dpu_bin_path DPU binary path
 * @param params NW parameters
 * @param ranks Number of ranks to use
 * @param set Dataset
 * @return std::vector<int>
 */
std::vector<int> dpu_16s_pipeline(std::filesystem::path dpu_bin_path, const NwParameters &p, size_t n_ranks, const Set &set);
#endif /* E6039E80_5D9F_462C_ACAE_D977B65797AC */
