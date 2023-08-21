/*
 * Copyright 2022 - UPMEM
 */

#ifndef EB929A76_D158_433E_98EA_D3AA35D04259
#define EB929A76_D158_433E_98EA_D3AA35D04259

#include "wram_aligned_buffer.h"

///////////////////////////////
//          64
///////////////////////////////

typedef struct mram_sequential_writer_64
{
    uint8_t *current;
    __mram_ptr uint8_t *mram_ptr;
} __attribute__((aligned(8))) mram_sequential_writer_64;

static inline mram_sequential_writer_64 get_mram_sequential_writer_64(WramAligned64 *mseq, __mram_ptr uint8_t *mram_ptr)
{
    mram_sequential_writer_64 writer = {mseq->buffer + (64LU * me()), mram_ptr};

    return writer;
}

static inline void mram_sequential_writer_64_set(mram_sequential_writer_64 mseq, uint8_t value)
{
    *(mseq.current) = value;
}

static inline void mram_sequential_writer_64_flush(mram_sequential_writer_64 mseq)
{
    mram_write((void *)((uintptr_t)mseq.current & 0xFFFFFFC0), mseq.mram_ptr, 64);
}

static inline void mram_sequential_writer_64_next(mram_sequential_writer_64 *mseq)
{
    __asm__ volatile(
        "add %[cur], %[cur], 1, nc6, 64f;"
        "add %[cur], %[cur], -64;"
        "sdma %[cur], %[ptr], 7;"
        "add %[ptr], %[ptr], 64;"
        "64:"
        : [cur] "+r"(mseq->current),
          [ptr] "+r"(mseq->mram_ptr)
        :);
}

#endif /* EB929A76_D158_433E_98EA_D3AA35D04259 */
