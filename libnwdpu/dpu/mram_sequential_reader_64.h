/*
 * Copyright 2022 - UPMEM
 */

#ifndef B620F271_F2CD_4B11_BAF3_6857E043B135
#define B620F271_F2CD_4B11_BAF3_6857E043B135

#include "wram_aligned_buffer.h"

///////////////////////////////
//          64
///////////////////////////////

struct MramSeqReader
{
    /// @brief Data for the sequential reader
    uint8_t *current;             /// pointer in the aligned buffer
    __mram_ptr uint8_t *mram_ptr; /// pointer to the current 64 block of MRAM in the buffer
} /// Hoping the compiler can see the double load one day !
__attribute__((aligned(8)));

typedef struct MramSeqReader MramSeqReader;

/**
 * @brief Get the mram sequential reader 64 object
 *
 * @param mseq pointer to wram aligned buffer group aware
 * @param mram_ptr pointer to MRAM start of sequence
 * @param id group id
 * @return MramSeqReader
 */
static inline MramSeqReader get_mram_sequential_reader_64(WramAligned64 *mseq, __mram_ptr uint8_t *mram_ptr, uint32_t id)
{
    MramSeqReader reader = {mseq->buffer + (64 * id) + 63, mram_ptr - 64};
    return reader;
}

/**
 * @brief get next value, efficient buffer overflow check for cache.
 *
 * @param mseq
 * @return uint8_t*
 */
static inline uint8_t *mram_sequential_reader_64_next(MramSeqReader *mseq)
{
    __asm__ volatile(
        "add %[cur], %[cur], 1, nc6, 64f;"
        "add %[ptr], %[ptr], 64;"
        "add %[cur], %[cur], -64;"
        "ldma %[cur], %[ptr], 7;"
        "64:"
        : [cur] "+r"(mseq->current),
          [ptr] "+r"(mseq->mram_ptr)
        :);

    return mseq->current;
}

#endif /* B620F271_F2CD_4B11_BAF3_6857E043B135 */
