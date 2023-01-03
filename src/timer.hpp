/*
 * Copyright 2022 - UPMEM
 */

#ifndef F9C8281F_F477_471C_B5BE_4ACE15AE8324
#define F9C8281F_F477_471C_B5BE_4ACE15AE8324

#include <stdint.h>
#include <time.h>

struct Time : public timespec
{
    operator double() { return static_cast<double>(tv_sec) + static_cast<double>(tv_nsec) * 1e-9; }
};

/**
 * @brief Timer to count both cpu and wall time. Used to see core usage.
 *
 */
class Timer
{
    timespec wall_time_start{};
    timespec cpu_time_start{};

public:
    Timer()
    {
        clock_gettime(CLOCK_MONOTONIC_RAW, &wall_time_start);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpu_time_start);
    }

    double wall() const
    {
        timespec wall_time_stop{};
        clock_gettime(CLOCK_MONOTONIC_RAW, &wall_time_stop);

        return Time{{wall_time_stop.tv_sec - wall_time_start.tv_sec, wall_time_stop.tv_nsec - wall_time_start.tv_nsec}};
    }

    double
    cpu() const
    {
        timespec cpu_time_stop{};
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpu_time_stop);

        return Time{{cpu_time_stop.tv_sec - cpu_time_start.tv_sec, cpu_time_stop.tv_nsec - cpu_time_start.tv_nsec}};
    }

    void reset()
    {
        clock_gettime(CLOCK_MONOTONIC_RAW, &wall_time_start);
        clock_gettime(CLOCK_PROCESS_CPUTIME_ID, &cpu_time_start);
    }

    void print() const
    {
        printf("times:\n"
               "  wall:  %.2f\n"
               "  cpu:   %.2f\n"
               "  ratio: %.1f\n",
               wall(), cpu(), cpu() / wall());
    }
};

#endif /* F9C8281F_F477_471C_B5BE_4ACE15AE8324 */
