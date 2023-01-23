/*
 * Copyright 2022 - UPMEM
 */

#ifndef B2D58B33_5E71_4A6D_A79C_CF6C1A63DD94
#define B2D58B33_5E71_4A6D_A79C_CF6C1A63DD94

#include "mram_sequential_reader_64.h"
#include "../../cdefs.h"

/**
 * @brief a sequential adn reader. It is group aware.
 *
 */
typedef struct dna_reader
{
    uint32_t chunck;                       /// values are 2bits, 4 values are read at a time and kept in chunck.
    uint32_t offset;                       /// offset of current value.
    mram_sequential_reader_64 mram_reader; /// mram reader to read MRAM efficiently.
} __attribute__((aligned(8))) dna_reader;

/**
 * @brief Get the adn reader object, initialize the structure
 *
 * @param buffer aligned buffer (group aware)
 * @param mram_ptr mram pointer to the sequence beginning
 * @param group reader group
 * @return dna_reader
 */
static inline dna_reader get_dna_reader(wram_aligned_buffer_64 *buffer, __mram_ptr uint8_t *mram_ptr, uint32_t group)
{
    dna_reader reader = {0, 8, get_mram_sequential_reader_64(buffer, mram_ptr, group)};

    return reader;
}

/**
 * @brief get next nucteotide, does not check bounds
 *
 * @param reader self explenatory
 * @return uint8_t next nucleotide
 */
static inline uint8_t dna_reader_next(dna_reader *reader)
{
    uint32_t offset = reader->offset + 2;
    uint32_t chunck = reader->chunck;

    if (offset > 6)
    {
        chunck = *mram_sequential_reader_64_next(&(reader->mram_reader));
        offset = 0;
    }
    reader->chunck = chunck;
    reader->offset = offset;

    return (chunck >> offset) & 3;
}

/**
 * @brief get next nucleotide in adn reader if i<j otherwise return c
 *
 * @param adn dna_reader
 * @param i represent index in adn
 * @param l represent size of adn
 * @param c default value if i > l
 * @return uint8_t
 */
static inline uint8_t next_nucleotide(dna_reader *adn, uint32_t i, uint32_t l, uint8_t c)
{
    if (i < l)
        return dna_reader_next(adn);

    return c;
}

#endif /* B2D58B33_5E71_4A6D_A79C_CF6C1A63DD94 */
