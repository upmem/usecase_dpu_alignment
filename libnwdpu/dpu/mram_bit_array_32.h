/*
 * Copyright 2022 - UPMEM
 */

#ifndef CD9DA9ED_324F_4111_B265_C9B05CE433CF
#define CD9DA9ED_324F_4111_B265_C9B05CE433CF

#include "wram_aligned_buffer.h"

typedef struct mram_bit_array_32
{
    __mram_ptr uint8_t *base_ptr;
    uint32_t current_offset;
    uint8_t *wram_buffer;
    bool written;
} __attribute__((aligned(8))) mram_bit_array_32;

static inline mram_bit_array_32 create_mram_bit_array_32(WramAligned32 *buffer, __mram_ptr uint8_t *mram_ptr, uint32_t id)
{
    mram_bit_array_32 array = {mram_ptr, UINT32_MAX, &buffer->buffer[32LU * id], false};
    return array;
}

static inline uint8_t mram_bit_array_32_get(mram_bit_array_32 *array, uint32_t i)
{
    uint32_t offset = i / 256;
    uint32_t bit_index = i % 256;
    uint32_t byte_index = bit_index / 8;
    uint32_t bit_offset = bit_index % 8;

    if (array->current_offset != offset)
    {
        if (array->written)
        {
            mram_write(array->wram_buffer, array->base_ptr + array->current_offset * 32LU, 32);
            array->written = false;
        }
        mram_read(array->base_ptr + offset * 32LU, array->wram_buffer, 32);
        array->current_offset = offset;
    }

    return 1 & (array->wram_buffer[byte_index] >> bit_offset);
}

static inline void mram_bit_array_32_set(mram_bit_array_32 *array, uint32_t i, uint8_t value)
{
    uint32_t offset = i / 256;
    uint32_t bit_index = i % 256;
    uint32_t byte_index = bit_index / 8;
    uint32_t bit_offset = bit_index % 8;

    value <<= bit_offset;
    uint8_t mask = ~(1 << bit_offset);

    if (array->current_offset != offset)
    {
        if (array->written)
            mram_write(array->wram_buffer, array->base_ptr + array->current_offset * 32LU, 32);

        mram_read(array->base_ptr + offset * 32LU, array->wram_buffer, 32);
        array->current_offset = offset;
    }

    uint8_t byte = array->wram_buffer[byte_index];
    byte = byte & mask;
    array->wram_buffer[byte_index] = byte | value;
    array->written = true;
}

#endif /* CD9DA9ED_324F_4111_B265_C9B05CE433CF */
