/*
 * Copyright 2022 - UPMEM
 */

#ifndef F9C8281F_F477_471C_B5BE_4ACE15AE8324
#define F9C8281F_F477_471C_B5BE_4ACE15AE8324

#include <stdint.h>
#include <time.h>

/**
 * @brief Wrapper around timespec for easy convertion
 *
 */
struct Time : public timespec
{
    /// @brief Initialize time with a specific clocks
    /// @param clock_id
    Time(clockid_t clock_id) : timespec()
    {
        clock_gettime(clock_id, this);
    }

    /// @brief Converts timespec to double
    operator double() const
    {
        return static_cast<double>(tv_sec) + static_cast<double>(tv_nsec) * 1e-9;
    }

    /// @brief Returns the difference between to times in floating point seconds
    /// @param t2
    /// @return
    double operator-(const Time &t2) const { return static_cast<double>(*this) - static_cast<double>(t2); }
};

/**
 * @brief Timer to count both CPU and Wall time. Used to see core usage.
 *
 */
class Timer
{
    Time wall_time_start{CLOCK_MONOTONIC_RAW};
    Time cpu_time_start{CLOCK_PROCESS_CPUTIME_ID};

public:
    /// @brief Returns the Wall time (Real world time)
    /// @return
    double Wall() const
    {
        return Time{CLOCK_MONOTONIC_RAW} - wall_time_start;
    }

    /// @brief Returns the CPU time (affected by multithreading)
    /// @return
    double CPU() const
    {
        return Time{CLOCK_PROCESS_CPUTIME_ID} - cpu_time_start;
    }

    /// @brief Prints timers and the ratio between the two
    /// @param prefix
    void Print(const std::string &prefix = "") const
    {
        printf(
            (prefix + "times:\n" +
             prefix + "  Wall:  %.2f\n" +
             prefix + "  CPU:   %.2f\n" +
             prefix + "  ratio: %.1f\n")
                .c_str(),
            Wall(), CPU(), CPU() / Wall());
    }
};

#endif /* F9C8281F_F477_471C_B5BE_4ACE15AE8324 */
