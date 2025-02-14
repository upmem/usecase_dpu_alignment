/*
 * Copyright 2022 - UPMEM
 */

#include <alloc.h>
#include <perfcounter.h>

#include "dna_reader.h"
#include "assert.h"
#include "nw_common.h"
#include "mram_2bits_array_64.h"
#include "mram_buffered_array_64.h"

#ifndef PERF_COUNT_TYPE
#define PERF_COUNT_TYPE COUNT_CYCLES
#endif

__host NwMetadataDPU metadata;
__mram NwSequenceMetadataMram sequence_metadata;

__host ComparisonMetadata meta_index;
__mram NwScoreOutput output;

__dma_aligned uint8_t buf_av[NR_GROUPS][W_MAX];
__dma_aligned uint8_t buf_bv[NR_GROUPS][W_MAX];

void align_initialisations()
{
  const uint32_t pool_id = group();

  align_data[pool_id].l1 = sequence_metadata.lengths[align_data[pool_id].s1];
  align_data[pool_id].av = buf_av[pool_id];
  align_data[pool_id].i = init_dna1(&sequences[sequence_metadata.indexes[align_data[pool_id].s1]]);

  align_data[pool_id].l2 = sequence_metadata.lengths[align_data[pool_id].s2];
  align_data[pool_id].bv = buf_bv[pool_id];
  align_data[pool_id].j = init_dna2(&sequences[sequence_metadata.indexes[align_data[pool_id].s2]]);

  init_pv();
  init_ppv();
  init_fv();
  init_ev();

  align_data[pool_id].dir = RIGHT;
  align_data[pool_id].prev_dir = RIGHT;
}

/**
 * @brief Align sequences from its group.
 *        Done with adaptive band as defined here:
 *        https://www.biorxiv.org/content/10.1101/130633v2
 *
 * @return Alignment score
 */
int align()
{
  // initialize all buffers and values

  const uint32_t pool_id = group();

  align_initialisations();

  int32_t down = 0;

  // Main DP loop, one iteration computes one frontwave
  for (uint32_t d = 1; d < align_data[pool_id].l1 + align_data[pool_id].l2; d++)
  {

    align_data[pool_id].prev_dir = align_data[pool_id].dir;
    align_data[pool_id].dir = next_direction(align_data[pool_id].pv, align_data[pool_id].i, align_data[pool_id].l1, align_data[pool_id].j, align_data[pool_id].l2);

    if (align_data[pool_id].dir == DOWN)
    {
      parallel_sr();
      down++;
      align_data[pool_id].uv = align_data[pool_id].pv;
      align_data[pool_id].lv = align_data[pool_id].pv - 1;
      shift_right_if_previous_direction_is_down(align_data[pool_id].prev_dir, align_data[pool_id].ppv);
    }
    else
    {
      parallel_sl();
      align_data[pool_id].lv = align_data[pool_id].pv;
      align_data[pool_id].uv = align_data[pool_id].pv + 1;
      shift_left_if_previous_direction_is_right(align_data[pool_id].prev_dir, align_data[pool_id].ppv);
    }
    wait_shift();

    compute_affine_score();
    // compute_affine_score_slow();

    int32_t *tmpv = align_data[pool_id].pv;
    align_data[pool_id].pv = align_data[pool_id].ppv;
    align_data[pool_id].ppv = tmpv;
  }

  return align_data[pool_id].pv[(W_MAX >> 1) + (down - align_data[pool_id].l2)];
}

extern uint64_t nw_perf_cnt;

uint32_t seq1_id = 0;
uint32_t seq2_id = 1;
uint32_t score_offset = 0;

MUTEX_INIT(seq_id_mutex);
MUTEX_INIT(score_mutex);
BARRIER_INIT(start_barrier, NR_TASKLETS);
BARRIER_INIT(end_barrier, NR_GROUPS);

void next_pair()
{
  score_offset++;
  seq2_id++;

  if (seq2_id == meta_index.size)
  {
    seq1_id++;
    seq2_id = seq1_id + 1;
  }
}

static inline uint32_t sum_integers(uint32_t i)
{
  return i * (i - 1) / 2;
}

static inline uint32_t triangular_index(uint32_t i, uint32_t j, uint32_t n)
{
  return sum_integers(n) - sum_integers(n - i) + j - i - 1;
}

int main()
{
  // 1) set a global ID (under mutex) to define for each tasklet
  //  - what is the next pair of sequences to look at.

  if (me() == 0)
  {
    mem_reset();
    score_offset = 0;
    seq1_id = meta_index.start_row;
    seq2_id = meta_index.start_col;

    perfcounter_config(PERF_COUNT_TYPE, true);
  }
  tasklet_params[me()].start = (me() % 4) * 32;
  barrier_wait(&start_barrier);

  wait_for_work();

  while (true)
  {
    mutex_lock(seq_id_mutex);
    uint32_t local_score_offset = score_offset;
    uint32_t seq1 = seq1_id;
    uint32_t seq2 = seq2_id;
    next_pair();

    mutex_unlock(seq_id_mutex);

    if (local_score_offset >= meta_index.count)
      break;

    // set parameter for the alignment group
    const uint32_t pool_id = group();
    align_data[pool_id].s1 = seq1;
    align_data[pool_id].s2 = seq2;
    align_data[pool_id].s_off = local_score_offset;

    int score = align();
    mutex_lock(score_mutex);
    output.scores[local_score_offset] = score;
    mutex_unlock(score_mutex);
  }

  barrier_wait(&end_barrier);

  if (me() == 0)
    output.perf_counter = perfcounter_get();

  return 0;
}
