/*
 * Copyright 2022 - UPMEM
 */

#ifndef EFB491BB_CE51_45FE_BB8B_8CD42179622B
#define EFB491BB_CE51_45FE_BB8B_8CD42179622B

#include <stdint.h>

#define DPU_MAX_NUMBER_OF_SEQUENCES 1024LU                // Max of total sequences in dpu
#define DPU_MAX_NUMBER_OF_SEQUENCES_MRAM 16384LU          // Max of total sequences in dpu
#define SCORE_METADATA_MAX_NUMBER_OF_SCORES 4096LU        // Max number of total alignment in dpu
#define SCORE_METADATA_MAX_NUMBER_OF_SCORES_MRAM 131072LU // Max number of total alignment in dpu
#define SCORE_METADATA_MAX_NUMBER_OF_SET 36LU             // Max number of set in a dpu
#define SCORE_MAX_SEQUENCES_TOTAL_SIZE 3840000LU          // 4MB of MRAM for sequences
#define METADATA_MAX_NUMBER_OF_SCORES 4096LU              // Max number of pair alignment
#define MAX_CIGAR_SIZE 32000000LU                         // 32MB of MRAM for cigars
#define DPU_MAX_SEQUENCE_SIZE 80000LU                     // Is use for direction bit array
#define W_MAX 128LU                                       // Width of anti-diagonal use in dpu

// typedef uint16_t value_t;

/**
 * @brief Represent the trace of previous cell score
 *
 */
typedef enum
{
    /// @brief  Score origin possible values
    DMISS = 0,  /// score is from mismatch
    DMATCH = 1, /// score is from match
    LEFT = 2,   /// score is from left gap
    UP = 3,     /// score is from upper gap
} TraceValue;

/**
 * @brief Represent the direction of the next band to compute
 *
 */
typedef enum
{
    /// @brief Direction possible values
    DOWN = 0, /// band goes down
    RIGHT = 1 /// band goes right
} Direction;

/**
 * @brief Structure for data exchange between host and DPU.
 * Contains index and length of sequences. Sequence buffer is
 * send separatly to keep the structure small enough otherwise
 *  stack errors occurs.
 *
 */
typedef struct NwMetadataDPU
{
    /// @brief Collection of all common data needed for N&W
    uint32_t indexes[DPU_MAX_NUMBER_OF_SEQUENCES];       /// index of Nth sequence in sequence buffer
    uint16_t lengths[DPU_MAX_NUMBER_OF_SEQUENCES];       /// length of Nth sequence
    uint8_t set_sizes[SCORE_METADATA_MAX_NUMBER_OF_SET]; /// Size of each set sent. Max set size is 255
    uint32_t number_of_sets;                             /// number of total set sent
    int32_t match;                                       /// match score
    int32_t mismatch;                                    /// mismatch score
    int32_t gap_opening;                                 /// gap opening score
    int32_t gap_extension;                               /// gap extension score
    // int8_t pad[4];                                       /// padding for mram compliance
} NwMetadataDPU;

typedef struct NwSequenceMetadataMram
{
    uint32_t indexes[DPU_MAX_NUMBER_OF_SEQUENCES_MRAM]; /// index of Nth sequence in sequence buffer
    uint16_t lengths[DPU_MAX_NUMBER_OF_SEQUENCES_MRAM]; /// length of Nth sequence
} NwSequenceMetadataMram;

/**
 * @brief Represents the needed parameters to compute upper triangular comparison matrix
 *
 */
typedef struct ComparisonMetadata
{
    /// @brief All against all needed values
    uint32_t start_row; /// Starting row for comparison matrix
    uint32_t start_col; /// Starting column for comparison matrix
    uint32_t count;     /// How many comparison to do
    uint32_t size;      /// Total number of sequence to compare (comparison matrix size)
} ComparisonMetadata;

/**
 * @brief Structure for data send back from DPU to host.
 * Contains the perfcounter, score of each pair alignment
 * and lenght of all cigars. CIGARs are sent separatly.
 *
 */
typedef struct NwCigarOutput
{
    /// @brief Relevant data
    uint64_t perf_counter;                           /// performance counter, cycle or instruction can be change on dpu code size.
    int32_t scores[METADATA_MAX_NUMBER_OF_SCORES];   /// score of pair alignment
    uint16_t lengths[METADATA_MAX_NUMBER_OF_SCORES]; /// length of CIGARs
} NwCigarOutput;

/**
 * @brief Structure for data send back from DPU to host.
 * Contains the perfcounter and score of each comparison.
 *
 */
typedef struct NwScoreOutput
{
    /// @brief Relevant data
    uint64_t perf_counter;                                    /// performance counter
    int32_t scores[SCORE_METADATA_MAX_NUMBER_OF_SCORES_MRAM]; /// score of comparison
} NwScoreOutput;

#endif /* EFB491BB_CE51_45FE_BB8B_8CD42179622B */
