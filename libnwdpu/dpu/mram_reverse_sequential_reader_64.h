/*
 * Copyright 2022 - UPMEM
 */

#ifndef A6B11990_1709_4186_A8AA_55948ADF43A7
#define A6B11990_1709_4186_A8AA_55948ADF43A7

#include "wram_aligned_buffer.h"

///////////////////////////////
//          64
///////////////////////////////

typedef struct mram_reverse_sequential_reader_64
{
    uint8_t *current;
    __mram_ptr uint8_t *mram_ptr;
} __attribute__((aligned(8))) mram_reverse_sequential_reader_64;

static inline mram_reverse_sequential_reader_64 get_mram_reverse_sequential_reader_64(WramAligned64 *buf, __mram_ptr uint8_t *mram_ptr)
{
    uintptr_t offset = (uintptr_t)mram_ptr & 0x3F;
    mram_ptr = (__mram_ptr uint8_t *)((uintptr_t)mram_ptr & 0xFFFFFFC0);
    mram_reverse_sequential_reader_64 reader = {buf->buffer + (64LU * me()), mram_ptr};

    if (offset != 0)
        mram_read(reader.mram_ptr, reader.current, offset);
    reader.current += offset;
    return reader;
}

static inline uint8_t *mram_reverse_sequential_reader_64_prev(mram_reverse_sequential_reader_64 *mseq)
{
    mseq->current--;

    if (((uintptr_t)mseq->current % 64) == 63)
    {
        mseq->current++;
        mseq->mram_ptr -= 64;
        mram_read(mseq->mram_ptr, mseq->current, 64);
        mseq->current += 63;
    }

    return mseq->current;
}

#endif /* A6B11990_1709_4186_A8AA_55948ADF43A7 */
