/*
 * Copyright 2022 - UPMEM
 */

#ifndef C82CA1A0_84CB_453A_8847_504E9690AAA8
#define C82CA1A0_84CB_453A_8847_504E9690AAA8

#include <defs.h>
#include <mram.h>
#include <stdbool.h>

// max number of tasklets
#ifndef NR_TASKLETS
#define NR_TASKLETS (-1)
#endif

// if no NR_GROUPS is given, make one buffer per tasklet
#ifndef NR_GROUPS
#define NR_GROUPS NR_TASKLETS
#endif

/**
 * @brief Create an aligned wram buffer of 32 bytes for all tasklets/groups
 *
 */
typedef struct wram_aligned_buffer_32
{
    uint8_t buffer[32 * NR_GROUPS];
} __attribute__((aligned(32))) wram_aligned_buffer_32;

/**
 * @brief Create an aligned wram buffer of 64 bytes for all tasklets/groups
 *
 */
typedef struct wram_aligned_buffer_64
{
    uint8_t buffer[64 * NR_GROUPS];
} __attribute__((aligned(64))) wram_aligned_buffer_64;

/**
 * @brief Create an aligned wram buffer of 128 bytes for all tasklets/groups
 *
 */
typedef struct wram_aligned_buffer_128
{
    uint8_t buffer[128 * NR_GROUPS];
} __attribute__((aligned(128))) wram_aligned_buffer_128;

/**
 * @brief Create an aligned wram buffer of 256 bytes for all tasklets/groups
 *
 */
typedef struct wram_aligned_buffer_256
{
    uint8_t buffer[256 * NR_GROUPS];
} __attribute__((aligned(256))) wram_aligned_buffer_256;

#endif /* C82CA1A0_84CB_453A_8847_504E9690AAA8 */
