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
#define NR_TASKLETS (24)
#endif

// if no NR_GROUPS is given, make one buffer per tasklet
#ifndef NR_GROUPS
#define NR_GROUPS NR_TASKLETS
#endif

/**
 * @brief Create an aligned wram buffer of 32 bytes for all tasklets/groups
 *
 */
struct WramAligned32
{
    /// @brief Aligned buffer Containing the data
    uint8_t buffer[32 * NR_GROUPS];
} /// Force alignment on 32 Bytes
__attribute__((aligned(32)));

/**
 * @brief Create an aligned wram buffer of 64 bytes for all tasklets/groups
 *
 */
struct WramAligned64
{
    /// @brief Aligned buffer Containing the data
    uint8_t buffer[64 * NR_GROUPS];
} /// Force alignment on 64 Bytes
__attribute__((aligned(64)));

/**
 * @brief Create an aligned wram buffer of 128 bytes for all tasklets/groups
 *
 */
struct WramAligned128
{
    /// @brief Aligned buffer Containing the data
    uint8_t buffer[128 * NR_GROUPS];
} /// Force alignment on 128 Bytes
__attribute__((aligned(128)));

/**
 * @brief Create an aligned wram buffer of 256 bytes for all tasklets/groups
 *
 */
struct WramAligned256
{
    /// @brief Aligned buffer Containing the data
    uint8_t buffer[256 * NR_GROUPS];
} /// Force alignment on 256 Bytes
__attribute__((aligned(256)));

/// @brief typedefs for easy use
typedef struct WramAligned32 WramAligned32;
typedef struct WramAligned64 WramAligned64;
typedef struct WramAligned128 WramAligned128;
typedef struct WramAligned256 WramAligned256;

#endif /* C82CA1A0_84CB_453A_8847_504E9690AAA8 */
