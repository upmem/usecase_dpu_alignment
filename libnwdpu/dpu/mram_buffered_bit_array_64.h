/*
 * Copyright 2022 - UPMEM
 */

#ifndef F1F34886_F528_48B3_AE83_B6A94E336779
#define F1F34886_F528_48B3_AE83_B6A94E336779

#include "wram_aligned_buffer.h"

typedef struct mram_buffered_bit_array_64
{
    __mram_ptr uint8_t *const base_ptr;
    uint32_t current_offset;
    uint8_t *const wram_buffer;
    bool written;
} __attribute__((aligned(8))) mram_buffered_bit_array_64;

static inline mram_buffered_bit_array_64 get_mram_buffered_bit_array_64(WramAligned64 *buffer, __mram_ptr uint8_t *mram_ptr)
{
    mram_buffered_bit_array_64 array = {mram_ptr, UINT32_MAX, (uint8_t *)buffer + (64LU * me()), false};
    return array;
}

static inline uint8_t mram_buffered_bit_array_64_get(mram_buffered_bit_array_64 *array, uint32_t i)
{
    uint32_t offset = i / 512;
    uint32_t bit_index = i % 512;
    uint32_t byte_index = bit_index / 8;
    uint32_t bit_offset = bit_index % 8;

    if (array->current_offset != offset)
    {
        if (array->written)
        {
            mram_write(array->wram_buffer, array->base_ptr + array->current_offset * 64, 64);
            array->written = false;
        }
        mram_read(array->base_ptr + offset * 64, array->wram_buffer, 64);
        array->current_offset = offset;
    }

    return 1 & (array->wram_buffer[byte_index] >> bit_offset);
}

static inline void mram_buffered_bit_array_64_set(mram_buffered_bit_array_64 *array, uint32_t i, uint8_t value)
{
    uint32_t offset = i / 512;
    uint32_t bit_index = i % 512;
    uint32_t byte_index = bit_index / 8;
    uint32_t bit_offset = bit_index % 8;

    value <<= bit_offset;
    uint8_t mask = ~(1 << bit_offset);

    if (array->current_offset != offset)
    {
        if (array->written)
            mram_write(array->wram_buffer, array->base_ptr + array->current_offset * 64, 64);

        mram_read(array->base_ptr + offset * 64, array->wram_buffer, 64);
        array->current_offset = offset;
    }

    uint8_t byte = array->wram_buffer[byte_index];
    byte = byte & mask;
    array->wram_buffer[byte_index] = byte | value;
    array->written = true;
}

#endif /* F1F34886_F528_48B3_AE83_B6A94E336779 */
