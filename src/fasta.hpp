/*
 * Copyright 2022 - UPMEM
 */

#ifndef D6AA45B0_8C2F_470E_9570_4BF160F4C5EE
#define D6AA45B0_8C2F_470E_9570_4BF160F4C5EE

#include "types.hpp"

/**
 * @brief Reads a fasta file with all sequences having a set id in comment
 *
 * @param filename
 * @return A vector of set
 */
Sets read_set_fasta(const std::filesystem::path &filename);

/**
 * @brief Reads all sequences in a fasta file
 *
 * @param filename
 * @return A vector of sequences
 */
Set read_seq_fasta(const std::filesystem::path &filename);

#endif /* D6AA45B0_8C2F_470E_9570_4BF160F4C5EE */
