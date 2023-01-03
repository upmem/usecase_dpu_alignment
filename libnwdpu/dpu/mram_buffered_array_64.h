/*
 * Copyright 2022 - UPMEM
 */

#ifndef CAB36ED5_72A0_49FC_89C9_9217059FAFCB
#define CAB36ED5_72A0_49FC_89C9_9217059FAFCB

#include "wram_aligned_buffer.h"

typedef struct mram_buffered_array_64
{
    __mram_ptr uint8_t *const base_ptr;
    uint32_t current_offset;
    uint8_t *const wram_buffer;
    bool written;
} __attribute__((aligned(8))) mram_buffered_array_64;

static inline mram_buffered_array_64 get_mram_buffered_array_64(wram_aligned_buffer_64 *buffer, __mram_ptr uint8_t *mram_ptr, uint32_t id)
{
    mram_buffered_array_64 array = {mram_ptr, UINT32_MAX, (uint8_t *)buffer + (64LU * id), false};
    return array;
}

static inline uint8_t mram_buffered_array_64_get(mram_buffered_array_64 *array, uint32_t i)
{
    uint32_t offset = i / 64;

    if (array->current_offset == offset)
        return array->wram_buffer[i % 64];

    if (array->written)
    {
        mram_write(array->wram_buffer, array->base_ptr + array->current_offset * 64, 64);
        array->written = false;
    }
    array->current_offset = offset;
    offset *= 64;

    mram_read(array->base_ptr + offset, array->wram_buffer, 64);

    return array->wram_buffer[i % 64];
}

static inline void mram_buffered_array_64_set(mram_buffered_array_64 *array, uint32_t i, uint8_t value)
{
    uint32_t offset = i / 64;

    if (array->current_offset != offset)
    {
        if (array->written)
            mram_write(array->wram_buffer, array->base_ptr + array->current_offset * 64, 64);

        mram_read(array->base_ptr + offset * 64, array->wram_buffer, 64);
        array->current_offset = offset;
    }

    array->wram_buffer[i % 64] = value;
    array->written = true;
}

static inline void mram_buffered_array_64_flush(mram_buffered_array_64 *array)
{
    if (array->written)
        mram_write(array->wram_buffer, array->base_ptr + array->current_offset * 64, 64);
}

#endif /* CAB36ED5_72A0_49FC_89C9_9217059FAFCB */
