/*
 * Copyright 2022 - UPMEM
 */

#ifndef A2FB4D02_DF1F_48BC_B3B1_0EA01AE7E160
#define A2FB4D02_DF1F_48BC_B3B1_0EA01AE7E160

#include "wram_aligned_buffer.h"

typedef struct mram_2bits_array_64
{
    __mram_ptr uint8_t *const base_ptr;
    uint32_t current_offset;
    uint8_t *const wram_buffer;
    bool written;
} __attribute__((aligned(8))) mram_2bits_array_64;

static inline mram_2bits_array_64 create_mram_2bits_array_64(WramAligned64 *buffer, __mram_ptr uint8_t *mram_ptr, uint32_t id)
{
    mram_2bits_array_64 array = {mram_ptr, UINT32_MAX, &buffer->buffer[64LU * id], false};
    return array;
}

static inline uint8_t mram_2bits_array_64_get(mram_2bits_array_64 *array, uint32_t i)
{
    uint32_t offset = i / 256;
    uint32_t bit_index = i % 256;
    uint32_t byte_index = bit_index / 4;
    uint32_t bit_offset = bit_index % 4;

    if (array->current_offset != offset)
    {
        if (array->written)
        {
            mram_write(array->wram_buffer, array->base_ptr + array->current_offset * 64LU, 64);
            array->written = false;
        }
        mram_read(array->base_ptr + offset * 64LU, array->wram_buffer, 64);
        array->current_offset = offset;
    }

    return 3 & (array->wram_buffer[byte_index] >> (bit_offset * 2));
}

static inline void mram_2bits_array_64_set(mram_2bits_array_64 *array, uint32_t i, uint8_t value)
{
    uint32_t offset = i / 256;
    uint32_t bit_index = i % 256;
    uint32_t byte_index = bit_index / 4;
    uint32_t bit_offset = bit_index % 4;

    bit_offset *= 2;
    value <<= bit_offset;
    uint8_t mask = ~(3 << bit_offset);

    if (array->current_offset != offset)
    {
        if (array->written)
            mram_write(array->wram_buffer, array->base_ptr + array->current_offset * 64LU, 64);

        mram_read(array->base_ptr + offset * 64LU, array->wram_buffer, 64);
        array->current_offset = offset;
    }

    uint8_t byte = array->wram_buffer[byte_index];
    byte = byte & mask;
    array->wram_buffer[byte_index] = byte | value;
    array->written = true;
}

#endif /* A2FB4D02_DF1F_48BC_B3B1_0EA01AE7E160 */
