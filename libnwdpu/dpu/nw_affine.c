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

__host NW_dpu_metadata_input metadata;
__mram uint8_t sequences[SCORE_MAX_SEQUENCES_TOTAL_SIZE];

__mram uint8_t dirs[NR_GROUPS][DPU_MAX_SEQUENCE_SIZE * 2 / 8]; // Bit array

__host NW_dpu_output output;
__mram uint8_t cigars[MAX_CIGAR_SIZE];
__mram uint32_t cigar_indexes[METADATA_MAX_NUMBER_OF_SCORES];

wram_aligned_buffer_64 dna_reader_buffer1;
wram_aligned_buffer_64 dna_reader_buffer2;

wram_aligned_buffer_32 direction_buffer;
wram_aligned_buffer_32 trace_wram_buffer;

wram_aligned_buffer_32 t_e_wram_buffer;
wram_aligned_buffer_32 t_f_wram_buffer;

__host struct m_buf
{
  __attribute__((aligned(64))) int32_t ev[W_MAX];
  __attribute__((aligned(64))) int32_t fv[W_MAX];
  __attribute__((aligned(64))) int32_t pv[W_MAX + 4];
  __attribute__((aligned(64))) int32_t ppv[W_MAX + 4];
} align_buffers[NR_GROUPS];

__dma_aligned uint8_t buf_av[NR_GROUPS][W_MAX];
__dma_aligned uint8_t buf_bv[NR_GROUPS][W_MAX];

void reverse(__mram_ptr uint8_t *mram, size_t size)
{
  for (size_t i = 0; i < size / 2; i++)
  {
    uint8_t tmp = mram[size - 1 - i];
    mram[size - 1 - i] = mram[i];
    mram[i] = tmp;
  }
}

void init_pv()
{
  const uint32_t align_id = group();
  const int32_t gapo = metadata.gap_opening;
  const int32_t gape = metadata.gap_extension;
  const int32_t gapoe = gapo + gape;

  align_data[align_id].pv = align_buffers[align_id].pv + 2;

  int32_t *pv = align_data[align_id].pv;
  pv[-1] = INT32_MIN / 2;
  pv[W_MAX] = INT32_MIN / 2;

  for (uint32_t i = 0; i < W_MAX; i++)
    pv[i] = INT32_MIN / 2;

  pv[W_MAX / 2] = -gapoe;
  pv[(W_MAX / 2) - 1] = -gapoe;
}

void init_ppv()
{
  const uint32_t align_id = group();
  align_data[align_id].ppv = align_buffers[align_id].ppv + 2;

  int32_t *ppv = align_data[align_id].ppv;
  ppv[-1] = INT32_MIN / 2;
  ppv[W_MAX] = INT32_MIN / 2;

  for (uint32_t i = 0; i < W_MAX; i++)
    ppv[i] = INT32_MIN / 2;

  ppv[W_MAX / 2] = 0;
}

void init_ev()
{
  const uint32_t align_id = group();
  const int32_t gapo = metadata.gap_opening;
  const int32_t gape = metadata.gap_extension;
  const int32_t gapoe = gapo + gape;

  align_data[align_id].ev = align_buffers[align_id].ev;

  int32_t *ev = align_data[align_id].ev;
  for (uint32_t i = 0; i < W_MAX; i++)
    ev[i] = INT32_MIN / 2;

  ev[(W_MAX / 2) - 1] = -gapoe;
}

void init_fv()
{
  const uint32_t align_id = group();
  const int32_t gapo = metadata.gap_opening;
  const int32_t gape = metadata.gap_extension;
  const int32_t gapoe = gapo + gape;

  align_data[align_id].fv = align_buffers[align_id].fv;

  int32_t *fv = align_data[align_id].fv;
  for (uint32_t i = 0; i < W_MAX; i++)
    fv[i] = INT32_MIN / 2;

  fv[(W_MAX / 2)] = -gapoe;
}

uint32_t init_dna1()
{
  const uint32_t align_id = group();

  align_data[align_id].dna1 = get_dna_reader(
      &dna_reader_buffer1,
      &sequences[metadata.indexes[align_data[align_id].s1]],
      align_id);

  int w2 = (W_MAX >> 1);

  dna_reader *s1 = &align_data[align_id].dna1;
  uint8_t *av = align_data[align_id].av;
  uint32_t l1 = align_data[align_id].l1;

  uint32_t i = 0;
  for (; i < w2; i++)
  {
    av[i] = 'X';
    if (i < l1)
      av[w2 + i] = dna_reader_next(s1);
    else
      av[w2 + i] = 'X';
  }
  return i;
}

uint32_t init_dna2()
{
  const uint32_t align_id = group();
  int w2 = (W_MAX >> 1);

  align_data[align_id].dna2 = get_dna_reader(
      &dna_reader_buffer2,
      &sequences[metadata.indexes[align_data[align_id].s2]],
      align_id);

  dna_reader *s2 = &align_data[align_id].dna2;
  uint8_t *bv = align_data[align_id].bv;
  uint32_t l2 = align_data[align_id].l2;

  uint32_t i = 0;
  for (; i < w2; i++)
  {
    bv[w2 + i] = 'Y';
    if (i < l2)
      bv[w2 - 1 - i] = dna_reader_next(s2);
    else
      bv[w2 - 1 - i] = 'Y';
  }
  return i;
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

  const uint32_t align_id = group();

  align_data[align_id].l1 = metadata.lengths[align_data[align_id].s1];
  align_data[align_id].av = buf_av[align_id];
  align_data[align_id].i = init_dna1();

  align_data[align_id].l2 = metadata.lengths[align_data[align_id].s2];
  align_data[align_id].bv = buf_bv[align_id];
  align_data[align_id].j = init_dna2();

  init_pv();
  init_ppv();
  init_fv();
  init_ev();

  align_data[align_id].dir = RIGHT;
  align_data[align_id].prev_dir = RIGHT;

  align_data[align_id].direction_array = create_mram_bit_array_32(&direction_buffer, dirs[align_id], align_id);

  align_data[align_id].down = 0;

  align_data[align_id].trace = trace_wram_buffer.buffer + (32LU * align_id);
  align_data[align_id].t_e = t_e_wram_buffer.buffer + (32LU * align_id);
  align_data[align_id].t_f = t_f_wram_buffer.buffer + (32LU * align_id);

  mram_bit_array_32_set(&align_data[align_id].direction_array, 0, RIGHT);

  align_data[align_id].trace[(W_MAX >> 1) / 4] = LEFT;
  align_data[align_id].trace[(W_MAX >> 1) / 4 - 1] = (UP << 6);
  mram_write(align_data[align_id].trace, trace_buffer[align_id], 32LU);

  for (int k = 0; k < W_MAX / 8; k++)
  {
    align_data[align_id].t_e[k] = 0xFF;
    align_data[align_id].t_f[k] = 0xFF;
  }

  mram_write(align_data[align_id].t_e, te_buffer[align_id], 16LU);
  mram_write(align_data[align_id].t_f, tf_buffer[align_id], 16LU);

  align_data[align_id].offset = W_MAX / 4;

  // Main DP loop, one iteration computes one frontwave
  for (uint32_t d = 1; d < align_data[align_id].l1 + align_data[align_id].l2; d++)
  {

    align_data[align_id].prev_dir = align_data[align_id].dir;
    align_data[align_id].dir = next_direction(align_data[align_id].pv, align_data[align_id].i, align_data[align_id].l1, align_data[align_id].j, align_data[align_id].l2);

    if (align_data[align_id].dir == DOWN)
    {
      parallel_sr();
      align_data[align_id].down++;
      align_data[align_id].uv = align_data[align_id].pv;
      align_data[align_id].lv = align_data[align_id].pv - 1;
      shift_right_if_previous_direction_is_down(align_data[align_id].prev_dir, align_data[align_id].ppv);
    }
    else
    {
      parallel_sl();
      align_data[align_id].lv = align_data[align_id].pv;
      align_data[align_id].uv = align_data[align_id].pv + 1;
      shift_left_if_previous_direction_is_right(align_data[align_id].prev_dir, align_data[align_id].ppv);
    }
    mram_bit_array_32_set(&align_data[align_id].direction_array, d, align_data[align_id].dir);
    wait_shift();

    compute_affine();

    mram_write(align_data[align_id].trace, trace_buffer[align_id] + align_data[align_id].offset, 32LU);
    mram_write(align_data[align_id].t_e, te_buffer[align_id] + (align_data[align_id].offset / 2), 16LU);
    mram_write(align_data[align_id].t_f, tf_buffer[align_id] + (align_data[align_id].offset / 2), 16LU);
    align_data[align_id].offset += 32LU;

    int32_t *tmpv = align_data[align_id].pv;
    align_data[align_id].pv = align_data[align_id].ppv;
    align_data[align_id].ppv = tmpv;
  }

  ///// Backtracing /////

  int32_t d = align_data[align_id].l1 + align_data[align_id].l2 - 1;

  uint32_t offset = ((align_data[align_id].l1 + align_data[align_id].l2) * W_MAX) - W_MAX + (W_MAX >> 1) + (align_data[align_id].down - align_data[align_id].l2);

  mram_buffered_array_64 res = get_mram_buffered_array_64(
      &dna_reader_buffer1,
      &cigars[cigar_indexes[align_data[align_id].s_off]],
      align_id);
  mram_2bits_array_64 trace_reader = create_mram_2bits_array_64(&dna_reader_buffer2, trace_buffer[align_id], align_id);
  mram_bit_array_32 te_reader = create_mram_bit_array_32(&t_e_wram_buffer, te_buffer[align_id], align_id);
  mram_bit_array_32 tf_reader = create_mram_bit_array_32(&t_f_wram_buffer, tf_buffer[align_id], align_id);

  uint32_t sp = 0; // number of steps, gives the cigar final size.
  for (; d >= 0; sp++)
  {

    Direction direction = mram_bit_array_32_get(&align_data[align_id].direction_array, d);
    Direction direction_prev;
    if (d > 0)
      direction_prev = mram_bit_array_32_get(&align_data[align_id].direction_array, d - 1);

    int o = (direction == RIGHT) ? 0 : 1;
    int o2 = (direction != direction_prev)
                 ? 0
                 : ((direction == RIGHT)
                        ? -1
                        : +1);
    uint8_t current_trace = mram_2bits_array_64_get(&trace_reader, offset);

    switch (current_trace)
    {
    case DMATCH:
      mram_buffered_array_64_set(&res, sp, '=');
      offset -= 2 * W_MAX + o2, d--;
      break;

    case DMISS:
      mram_buffered_array_64_set(&res, sp, 'X');
      offset -= 2 * W_MAX + o2, d--;
      break;

    case LEFT:
      // if gap, need to go back up to the gap beginning.
      while (mram_bit_array_32_get(&tf_reader, offset) == 0)
      {
        mram_buffered_array_64_set(&res, sp, 'I');
        offset -= W_MAX + o;
        d--;
        sp++;
        if (d < 0)
          break;
        direction = mram_bit_array_32_get(&align_data[align_id].direction_array, d);
        o = (direction == RIGHT) ? 0 : 1;
      }

      mram_buffered_array_64_set(&res, sp, 'I');
      offset -= W_MAX + o;
      break;

    case UP:
      // if gap, need to go back up to the gap beginning.
      while (mram_bit_array_32_get(&te_reader, offset) == 0)
      {
        mram_buffered_array_64_set(&res, sp, 'D');
        offset -= W_MAX - 1 + o;
        d--;
        sp++;
        if (d < 0)
          break;
        direction = mram_bit_array_32_get(&align_data[align_id].direction_array, d);
        o = (direction == RIGHT) ? 0 : 1;
      }
      mram_buffered_array_64_set(&res, sp, 'D');
      offset -= W_MAX - 1 + o;
      break;
    }

    d--;
  }

  output.lengths[align_data[align_id].s_off] = sp;

  mram_buffered_array_64_flush(&res);

  // traceback is from end to start. cigar needs to be change to start to end.
  reverse(&cigars[cigar_indexes[align_data[align_id].s_off]], sp);

  return align_data[align_id].pv[(W_MAX >> 1) + (align_data[align_id].down - align_data[align_id].l2)];
}

extern uint64_t nw_perf_cnt;

uint32_t set_id = 0;
uint32_t seq1_id = 0;
uint32_t seq2_id = 1;
uint32_t score_offset = 0;
uint32_t set_offset = 0;

MUTEX_INIT(seq_id_mutex);
BARRIER_INIT(start_barrier, NR_TASKLETS);
BARRIER_INIT(end_barrier, NR_GROUPS);

void next_pair()
{
  score_offset++;
  seq2_id++;

  if (seq2_id == metadata.set_sizes[set_id])
  {
    seq1_id++;
    seq2_id = seq1_id + 1;
  }
  if (seq1_id == metadata.set_sizes[set_id] - 1)
  {
    seq1_id = 0;
    seq2_id = 1;
    set_offset += metadata.set_sizes[set_id];
    set_id++;
  }
}

int main()
{
  // 1) set a global ID (under mutex) to define for each tasklet
  //  - what is the next pair of sequences to look at.

  if (me() == 0)
  {
    mem_reset();
    set_id = 0;
    score_offset = 0;
    set_offset = 0;
    seq1_id = 0;
    seq2_id = 1;

    for (int id = 0; id < NR_TASKLETS; id++)
      tasklet_params[id].start = (id % 4) * 32;
    perfcounter_config(PERF_COUNT_TYPE, true);
  }
  barrier_wait(&start_barrier);

  wait_for_work();

  while (true)
  {
    mutex_lock(seq_id_mutex);
    uint32_t local_score_offset = score_offset;
    uint32_t local_set_offset = set_offset;
    uint32_t local_set_id = set_id;
    uint32_t seq1 = seq1_id;
    uint32_t seq2 = seq2_id;

    next_pair();

    mutex_unlock(seq_id_mutex);

    if (local_set_id >= metadata.number_of_sets)
      break;

    // set parameter for the alignment group
    const uint32_t align_id = group();
    align_data[align_id].s1 = local_set_offset + seq1;
    align_data[align_id].s2 = local_set_offset + seq2;
    align_data[align_id].s_off = local_score_offset;

    output.scores[local_score_offset] = align();
  }

  barrier_wait(&end_barrier);

  if (me() == 0)
    output.perf_counter = perfcounter_get();

  return 0;
}
