/*
 * Copyright 2022 - UPMEM
 */

#ifndef AC0C563D_AFD5_4A05_9BF9_F00902ACD05C
#define AC0C563D_AFD5_4A05_9BF9_F00902ACD05C

#include <built_ins.h>
#include <mutex.h>
#include <profiling.h>
#include <stddef.h>

#include <barrier.h>

#include "mram_bit_array_32.h"

/**
 * @brief Structure of all variable needed to be shared in a group.
 * Naming follow the paper: https://www.biorxiv.org/content/early/2017/09/07/130633.full.pdf
 *
 */
struct align_t
{
    uint32_t s1;                       /// first sequence index
    uint32_t s2;                       /// second sequence index
    uint32_t l1;                       /// first sequence lenght
    uint32_t l2;                       /// second sequence length
    dna_reader dna1;                   /// first sequence dna reader
    dna_reader dna2;                   /// second sequence dna reader
    uint32_t s_off;                    /// index of results: score and cigar
    uint8_t *av;                       /// pointer to av buffer 128 value
    uint8_t *bv;                       /// pointer to bv buffer 128 value
    uint32_t i;                        /// current position on sequence 1
    uint32_t j;                        /// current position on sequence 2
    int32_t *pv;                       /// pointer to pv buffer 128 values. -1 and 128 are valid.
    int32_t *ppv;                      /// pointer to ppv buffer 128 values. -1 and 128 are valid.
    int32_t *ev;                       /// pointer to ev buffer
    int32_t *fv;                       /// pointer to fv buffer
    int32_t *uv;                       /// pointer to pv or pv+1
    int32_t *lv;                       /// pointer to pv-1 or pv
    uint8_t *trace;                    /// pointer to trace buffer
    Direction dir;                     /// current band direction
    Direction prev_dir;                /// previous band direction
    mram_bit_array_32 direction_array; /// bit array keeping all band direction in MRAM
    uint8_t *t_e;                      /// pointer to E trace, for gap extension during backtrace
    uint8_t *t_f;                      /// pointer to F trace, for gap extension during backtrace
} align_data[NR_GROUPS];

/**
 * @brief There is 24 tasklets split into 6 group. Each 4 consecutive tasklets are assign the same group.
 *
 * @return group of a tasklet
 */
static inline sysname_t group()
{
    return me() / 4;
}

__mram_noinit uint8_t trace_buffer[NR_GROUPS][40000 / 2 * W_MAX]; /// 2 * 40000 seq on 2bits
__mram_noinit uint8_t te_buffer[NR_GROUPS][40000 / 4 * W_MAX];    /// 2 * 40000 seq on 1bit
__mram_noinit uint8_t tf_buffer[NR_GROUPS][40000 / 4 * W_MAX];    /// 2 * 40000 seq on 1bit
__mram_noinit uint8_t sequences[SCORE_MAX_SEQUENCES_TOTAL_SIZE];
WramAligned64 dna_reader_buffer1;
WramAligned64 dna_reader_buffer2;

extern NwMetadataDPU metadata;

/**
 * @brief shift right a 128 bytes buffer
 *
 * @param vec
 */
static inline void shift_right_u8(uint8_t *vec)
{
    uint64_t *p = (uint64_t *)vec;

#pragma unroll
    for (int i = 15; i > 0; i--)
    {
        p[i] = p[i] << 8;
        ((uint8_t *)(p + i))[0] = ((uint8_t *)(p + i))[-1];
    }
    p[0] = p[0] << 8;
}

/**
 * @brief shift left a 128 bytes buffer
 *
 * @param vec
 */
static inline void shift_left_u8(uint8_t *vec)
{
    uint64_t *p = (uint64_t *)vec;

#pragma unroll
    for (uint32_t i = 0; i < 15; i++)
    {
        p[i] = p[i] >> 8;
        ((uint8_t *)(p + i))[7] = ((uint8_t *)(p + i))[8];
    }
    p[15] = p[15] >> 8;
}

/**
 * @brief shift right a 128 int32_t buffer
 *
 * @param vec
 */
static inline void shift_right_s(int32_t *vec)
{
    uint64_t *p = (uint64_t *)vec;

    vec[127] = vec[126];

#pragma unroll
    for (int i = 62; i >= 0; i--)
    {
        register uint64_t val = p[i];
        vec[i * 2 + 1] = (int32_t)(val);
        vec[i * 2 + 2] = (int32_t)(val >> 32);
    }
}

/**
 * @brief shift left a 128 int32_t buffer
 *
 * @param vec
 */
static inline void shift_left_s(int32_t *vec)
{
    uint64_t *p = (uint64_t *)vec;

    vec[0] = vec[1];
#pragma unroll
    for (int i = 1; i < 64; i++)
    {
        register uint64_t val = p[i];
        vec[i * 2] = (int32_t)(val >> 32);
        vec[i * 2 - 1] = (int32_t)val;
    }
}

__host struct m_buf
{
    __attribute__((aligned(64))) int32_t ev[W_MAX];
    __attribute__((aligned(64))) int32_t fv[W_MAX];
    __attribute__((aligned(64))) int32_t pv[W_MAX + 4];
    __attribute__((aligned(64))) int32_t ppv[W_MAX + 4];
} align_buffers[NR_GROUPS];

static void init_pv()
{
    const uint32_t pool_id = group();
    const int32_t gapo = metadata.gap_opening;
    const int32_t gape = metadata.gap_extension;
    const int32_t gapoe = gapo + gape;

    align_data[pool_id].pv = align_buffers[pool_id].pv + 2;

    int32_t *pv = align_data[pool_id].pv;
    pv[-1] = INT32_MIN / 2;
    pv[W_MAX] = INT32_MIN / 2;

    for (uint32_t i = 0; i < W_MAX; i++)
        pv[i] = INT32_MIN / 2;

    pv[W_MAX / 2] = -gapoe;
    pv[(W_MAX / 2) - 1] = -gapoe;
}

static void init_ppv()
{
    const uint32_t pool_id = group();
    align_data[pool_id].ppv = align_buffers[pool_id].ppv + 2;

    int32_t *ppv = align_data[pool_id].ppv;
    ppv[-1] = INT32_MIN / 2;
    ppv[W_MAX] = INT32_MIN / 2;

    for (uint32_t i = 0; i < W_MAX; i++)
        ppv[i] = INT32_MIN / 2;

    ppv[W_MAX / 2] = 0;
}

static void init_ev()
{
    const uint32_t pool_id = group();
    const int32_t gapo = metadata.gap_opening;
    const int32_t gape = metadata.gap_extension;
    const int32_t gapoe = gapo + gape;

    align_data[pool_id].ev = align_buffers[pool_id].ev;

    int32_t *ev = align_data[pool_id].ev;
    for (uint32_t i = 0; i < W_MAX; i++)
        ev[i] = INT32_MIN / 2;

    ev[(W_MAX / 2) - 1] = -gapoe;
}

static void init_fv()
{
    const uint32_t pool_id = group();
    const int32_t gapo = metadata.gap_opening;
    const int32_t gape = metadata.gap_extension;
    const int32_t gapoe = gapo + gape;

    align_data[pool_id].fv = align_buffers[pool_id].fv;

    int32_t *fv = align_data[pool_id].fv;
    for (uint32_t i = 0; i < W_MAX; i++)
        fv[i] = INT32_MIN / 2;

    fv[(W_MAX / 2)] = -gapoe;
}

static uint32_t init_dna1(__mram_ptr uint8_t *index)
{
    const uint32_t pool_id = group();

    align_data[pool_id].dna1 = get_dna_reader(
        &dna_reader_buffer1,
        index,
        pool_id);

    int w2 = (W_MAX >> 1);

    dna_reader *reader = &align_data[pool_id].dna1;
    uint8_t *av = align_data[pool_id].av;
    uint32_t l1 = align_data[pool_id].l1;

    uint32_t i = 0;
    for (; i < w2; i++)
    {
        av[i] = 'X';
        if (i < l1)
            av[w2 + i] = dna_reader_next(reader);
        else
            av[w2 + i] = 'X';
    }
    return i;
}

static uint32_t init_dna2(__mram_ptr uint8_t *index)
{
    const uint32_t pool_id = group();
    int w2 = (W_MAX >> 1);

    align_data[pool_id].dna2 = get_dna_reader(
        &dna_reader_buffer2,
        index,
        pool_id);

    dna_reader *reader = &align_data[pool_id].dna2;
    uint8_t *bv = align_data[pool_id].bv;
    uint32_t l2 = align_data[pool_id].l2;

    uint32_t i = 0;
    for (; i < w2; i++)
    {
        bv[w2 + i] = 'Y';
        if (i < l2)
            bv[w2 - 1 - i] = dna_reader_next(reader);
        else
            bv[w2 - 1 - i] = 'Y';
    }
    return i;
}

/**
 * @brief wrapper around the cmpb4 instruction
 *
 * @param cmp 4 bytes result
 * @param av pointer to vector of 4 bytes to compare, needs to be align on 32bits
 * @param bv pointer to vector of 4 bytes to compare, needs to be align on 32bits
 */
static void compare_4bytes(uint32_t *cmp, const uint8_t *av, const uint8_t *bv)
{
    __asm__("cmpb4 %[cmp], %[av], %[bv]"
            : [cmp] "=r"(*cmp)
            : [av] "r"(*(uint32_t *)av), [bv] "r"(*(uint32_t *)bv));
}

/**
 * @brief shift ppv right if previous direction was down
 *
 * @param prev_dir
 * @param ppv
 */
static inline void shift_right_if_previous_direction_is_down(Direction prev_dir, int32_t *ppv)
{
    if (prev_dir == DOWN)
        shift_right_s(ppv), ppv[0] = INT32_MIN / 2;
}

/**
 * @brief shift ppv left if previous direction was right
 *
 * @param prev_dir
 * @param ppv
 */
static inline void shift_left_if_previous_direction_is_right(Direction prev_dir, int32_t *ppv)
{
    if (prev_dir == RIGHT)
        shift_left_s(ppv), ppv[W_MAX - 1] = INT32_MIN / 2;
}

/**
 * @brief Gives you direction with the greatest current score
 *
 * @param pv current buffer
 * @param w width
 * @param i position in query
 * @param l1 query size
 * @param j position in target
 * @param l2 target size
 * @return Direction
 */
static inline Direction next_direction(const int32_t *pv, uint32_t i, uint32_t l1, uint32_t j, uint32_t l2)
{
    if ((pv[0] > pv[W_MAX - 1] || i >= l1) && j < l2)
        return DOWN;

    return RIGHT;
}

static inline void send_work();
static inline void send_work_score();
static inline void wait_for_work();
static inline void wait_empty_slaves();

enum Functions
{
    AFFINE,
    SCORE_AFFINE,
    SHIFT_BV,
    SHIFT_AV,
    SHIFT_EV,
    SHIFT_FV
};

/**
 * @brief Structure to send parameters to sleeping tasklets in group pool.
 *
 */
struct m_param
{
    uint32_t start;      /// starting point for score computation in compute_affine
    enum Functions func; /// function the tasklet needs to execute upon wake up
} tasklet_params[NR_TASKLETS];

/**
 * @brief Compute all values of the new band.
 *
 */
static inline void compute_affine()
{
    send_work();

    const uint32_t pool_id = group();

    // creating alias for readability, does not impact performance
    const uint8_t *av = align_data[pool_id].av;
    const uint8_t *bv = align_data[pool_id].bv;
    int32_t *ppv = align_data[pool_id].ppv;
    int32_t *lv = align_data[pool_id].lv;
    int32_t *uv = align_data[pool_id].uv;
    int32_t *ev = align_data[pool_id].ev;
    int32_t *fv = align_data[pool_id].fv;
    uint8_t *traces = align_data[pool_id].trace;
    uint8_t *t_e = align_data[pool_id].t_e;
    uint8_t *t_f = align_data[pool_id].t_f;

    int32_t gape = metadata.gap_extension;
    int32_t gapoe = metadata.gap_opening + gape;
    int32_t match = metadata.match;
    int32_t miss = metadata.mismatch;

    uint32_t tef = 0; // E and F traces on the same register, can shift both in one instruction
    int count = 0;    // No need to optimize it, compiler do it fine.

    for (uint32_t wi = tasklet_params[me()].start; wi < tasklet_params[me()].start + 32; wi += 4)
    {
        uint8_t trace = 0;
        uint32_t cmp;
        compare_4bytes(&cmp, &av[wi], &bv[wi]);
        // #pragma unroll
        for (int i = 0; i < 4; i++)
        {
            uint8_t t = 0;
            int32_t tmppv = ppv[wi + i];
            int32_t tmpev = ev[wi + i];
            int32_t tmpfv = fv[wi + i];
            int32_t tmpuv = uv[wi + i];
            int32_t tmplv = lv[wi + i];

            __asm__(
                "lsr %[cmp], %[cmp], 8, so, 1f;"
                "move %[t], 0;"
                "add %[tmppv], %[tmppv], %[miss], true, 2f;"
                "1:"
                "move %[t], 64;"
                "add %[tmppv], %[tmppv], %[match];"
                "2:"
                "sub %[tmpuv], %[tmpuv], %[gapoe];" // 1
                "sub %[tmpev], %[tmpev], %[gape];"  // 0
                "lsr %[tef], %[tef], 1;"
                "jgts %[tmpev], %[tmpuv], 3f;"
                "move %[tmpev], %[tmpuv];"
                "or %[tef], %[tef], 128;"
                "3:"
                "sub %[tmplv], %[tmplv], %[gapoe];"
                "sub %[tmpfv], %[tmpfv], %[gape];"
                "jgts %[tmpfv], %[tmplv], 4f;"
                "move %[tmpfv], %[tmplv];"
                "or %[tef], %[tef], 8388608;"
                "4:"
                "jges %[tmppv], %[tmpev], 5f;"
                "move %[tmppv], %[tmpev];"
                "move %[t], 192;"
                "5:"
                "jges %[tmppv], %[tmpfv], 6f;"
                "move %[tmppv], %[tmpfv];"
                "move %[t], 128;"
                "6:"
                : [t] "=r"(t),
                  [cmp] "+r"(cmp),
                  [tmppv] "+r"(tmppv),
                  [tmpev] "+r"(tmpev),
                  [tmpfv] "+r"(tmpfv),
                  [tmpuv] "+r"(tmpuv),
                  [tmplv] "+r"(tmplv),
                  [tef] "+r"(tef)
                : [match] "r"(match),
                  [miss] "r"(miss),
                  [gapoe] "r"(gapoe),
                  [gape] "r"(gape)
                :);

            trace = (trace >> 2) | t;
            ppv[wi + i] = tmppv;
            ev[wi + i] = tmpev;
            fv[wi + i] = tmpfv;
        }

        traces[wi / 4] = trace;

        if (count++ == 1)
        {
            t_e[wi / 8] = tef;
            t_f[wi / 8] = tef >> 16;
            tef = 0;
            count = 0;
        }
    }

    wait_empty_slaves();
}

/**
 * @brief Compute all values of the new band.
 *
 */
static inline void compute_affine_score()
{
    send_work_score();

    const uint32_t pool_id = group();

    // creating alias for readability, does not impact performance
    const uint8_t *av = align_data[pool_id].av;
    const uint8_t *bv = align_data[pool_id].bv;
    int32_t *ppv = align_data[pool_id].ppv;
    int32_t *lv = align_data[pool_id].lv;
    int32_t *uv = align_data[pool_id].uv;
    int32_t *ev = align_data[pool_id].ev;
    int32_t *fv = align_data[pool_id].fv;

    int32_t gape = metadata.gap_extension;
    int32_t gapoe = metadata.gap_opening + gape;
    int32_t match = metadata.match;
    int32_t miss = metadata.mismatch;

    for (uint32_t wi = tasklet_params[me()].start; wi < tasklet_params[me()].start + 32; wi += 4)
    {
        uint32_t cmp;
        compare_4bytes(&cmp, &av[wi], &bv[wi]);
        // #pragma unroll
        for (int i = 0; i < 4; i++)
        {
            int32_t tmppv = ppv[wi + i];
            int32_t tmpev = ev[wi + i];
            int32_t tmpfv = fv[wi + i];
            int32_t tmpuv = uv[wi + i];
            int32_t tmplv = lv[wi + i];

            __asm__(
                "lsr %[cmp], %[cmp], 8, so, 1f;"
                "add %[tmppv], %[tmppv], %[miss], true, 2f;"
                "1:"
                "add %[tmppv], %[tmppv], %[match];"
                "2:"
                "sub %[tmpuv], %[tmpuv], %[gapoe];" // 1
                "sub %[tmpev], %[tmpev], %[gape];"  // 0
                "jgts %[tmpev], %[tmpuv], 3f;"
                "move %[tmpev], %[tmpuv];"
                "3:"
                "sub %[tmplv], %[tmplv], %[gapoe];"
                "sub %[tmpfv], %[tmpfv], %[gape];"
                "jgts %[tmpfv], %[tmplv], 4f;"
                "move %[tmpfv], %[tmplv];"
                "4:"
                "jges %[tmppv], %[tmpev], 5f;"
                "move %[tmppv], %[tmpev];"
                "5:"
                "jges %[tmppv], %[tmpfv], 6f;"
                "move %[tmppv], %[tmpfv];"
                "6:"
                : [cmp] "+r"(cmp),
                  [tmppv] "+r"(tmppv),
                  [tmpev] "+r"(tmpev),
                  [tmpfv] "+r"(tmpfv),
                  [tmpuv] "+r"(tmpuv),
                  [tmplv] "+r"(tmplv)
                : [match] "r"(match),
                  [miss] "r"(miss),
                  [gapoe] "r"(gapoe),
                  [gape] "r"(gape)
                :);

            ppv[wi + i] = tmppv;
            ev[wi + i] = tmpev;
            fv[wi + i] = tmpfv;
        }
    }

    wait_empty_slaves();
}

BARRIER_INIT(b_shift1, 3);
BARRIER_INIT(b_shift2, 3);
BARRIER_INIT(b_shift3, 3);
BARRIER_INIT(b_shift4, 3);
BARRIER_INIT(b_shift5, 3);
BARRIER_INIT(b_shift6, 3);

/**
 * @brief barrier to synchronize tasklets in a group after shifting buffers
 *
 */
static inline void wait_shift()
{
    if (group() == 0)
        barrier_wait(&b_shift1);
    else if (group() == 1)
        barrier_wait(&b_shift2);
    else if (group() == 2)
        barrier_wait(&b_shift3);
    else if (group() == 3)
        barrier_wait(&b_shift4);
    else if (group() == 4)
        barrier_wait(&b_shift5);
    else if (group() == 5)
        barrier_wait(&b_shift6);
}

/**
 * @brief shift the av buffer left then add next nucleotide
 *
 */
static inline void shift_av()
{
    const uint32_t pool_id = group();
    uint8_t *av = align_data[pool_id].av;

    shift_left_u8(av);
    av[W_MAX - 1] = next_nucleotide(
        &align_data[pool_id].dna1,
        align_data[pool_id].i++,
        align_data[pool_id].l1,
        'X');
}

/**
 * @brief shift bv buffer right and add next nucleotide
 *
 */
static inline void shift_bv()
{
    const uint32_t pool_id = group();
    uint8_t *bv = align_data[pool_id].bv;

    shift_right_u8(bv);
    bv[0] = next_nucleotide(
        &align_data[pool_id].dna2,
        align_data[pool_id].j++,
        align_data[pool_id].l2,
        'Y');
}

/**
 * @brief tasklets are waiting for next function to compute.
 * first tasklet of group exit early to become main.
 *
 */
static inline void wait_for_work()
{
    if (me() % 4 == 0)
        return;

    sysname_t master = group();

    while (true)
    {
        __asm__ volatile("stop;");

        switch (tasklet_params[me()].func)
        {
        case AFFINE:
            compute_affine();
            break;

        case SCORE_AFFINE:
            compute_affine_score();
            break;

        case SHIFT_BV:
            shift_bv();
            wait_shift();
            break;

        case SHIFT_FV:
            shift_right_s(align_data[master].fv);
            align_data[master].fv[0] = INT32_MIN / 2;
            wait_shift();
            break;

        case SHIFT_AV:
            shift_av();
            wait_shift();
            break;

        case SHIFT_EV:
            shift_left_s(align_data[master].ev);
            align_data[master].ev[W_MAX - 1] = INT32_MIN / 2;
            wait_shift();
            break;
        }
    }
}

/**
 * @brief resume tasklets, shift bv and fv buffers.
 *
 */
static inline void parallel_sr()
{
    sysname_t id1 = me() + 3;
    sysname_t id2 = me() + 2;

    tasklet_params[id1].func = SHIFT_BV;
    tasklet_params[id2].func = SHIFT_FV;

    __asm__ volatile("resume %[id], 0;" ::[id] "r"(id1));
    __asm__ volatile("resume %[id], 0;" ::[id] "r"(id2));
}

/**
 * @brief resume tasklets, shift av and ev buffers.
 *
 */
static inline void parallel_sl()
{
    sysname_t id1 = me() + 3;
    sysname_t id2 = me() + 2;

    tasklet_params[id1].func = SHIFT_AV;
    tasklet_params[id2].func = SHIFT_EV;

    __asm__ volatile("resume %[id], 0;" ::[id] "r"(id1));
    __asm__ volatile("resume %[id], 0;" ::[id] "r"(id2));
}

/**
 * @brief Wake up all group tasklets, compute the 4 part of new band.
 * Computation are independant.
 *
 */
static inline void send_work()
{
    if ((me() % 4) != 0)
        return;

    sysname_t id = me() + 1;
    tasklet_params[id].func = AFFINE;
    __asm__ volatile("resume %[id], 0;" ::[id] "r"(id));
    id++;
    tasklet_params[id].func = AFFINE;
    __asm__ volatile("resume %[id], 0;" ::[id] "r"(id));
    id++;
    tasklet_params[id].func = AFFINE;
    __asm__ volatile("resume %[id], 0;" ::[id] "r"(id));
}

static inline void send_work_score()
{
    if ((me() % 4) != 0)
        return;

    sysname_t id = me() + 1;
    tasklet_params[id].func = SCORE_AFFINE;
    __asm__ volatile("resume %[id], 0;" ::[id] "r"(id));
    id++;
    tasklet_params[id].func = SCORE_AFFINE;
    __asm__ volatile("resume %[id], 0;" ::[id] "r"(id));
    id++;
    tasklet_params[id].func = SCORE_AFFINE;
    __asm__ volatile("resume %[id], 0;" ::[id] "r"(id));
}

BARRIER_INIT(barrier1, 4);
BARRIER_INIT(barrier2, 4);
BARRIER_INIT(barrier3, 4);
BARRIER_INIT(barrier4, 4);
BARRIER_INIT(barrier5, 4);
BARRIER_INIT(barrier6, 4);

/**
 * @brief wait for all compute_affine tasklets to be finished for a group.
 *
 */
static inline void wait_empty_slaves()
{
    if (group() == 0)
        barrier_wait(&barrier1);
    else if (group() == 1)
        barrier_wait(&barrier2);
    else if (group() == 2)
        barrier_wait(&barrier3);
    else if (group() == 3)
        barrier_wait(&barrier4);
    else if (group() == 4)
        barrier_wait(&barrier5);
    else if (group() == 5)
        barrier_wait(&barrier6);
}

#endif /* AC0C563D_AFD5_4A05_9BF9_F00902ACD05C */
