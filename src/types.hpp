/*
 * Copyright 2022 - UPMEM
 */

#ifndef AD18B383_F97E_4512_98DA_46CE2947ACDD
#define AD18B383_F97E_4512_98DA_46CE2947ACDD

#include <array>
#include <filesystem>
#include <vector>

#include "../cdefs.h"
#include "timer.hpp"

/**
 * @brief Needleman & Wunsch parameters
 *
 */
struct NW_Parameters
{
    /// @brief Public parameters
    int match;             /// match bonus
    int mismatch;          /// mismatch penalty
    int32_t gap_opening;   /// gap opening penalty
    int32_t gap_extension; /// gap extension penalty
    int width;             /// band width

    /**
     * @brief Print Needleman & Wunsch parameters
     *
     */
    void print() const
    {
        printf("Alignment parameters:\n"
               "  match:         %d\n"
               "  mismatch:      %d\n"
               "  gap opening:   %d\n"
               "  gap extension: %d\n"
               "  width:         %d\n\n",
               match, mismatch, gap_opening, gap_extension, width);
    }
};

struct Cigar : public std::string
{
    int count_score(const NW_Parameters &params) const
    {
        int score = 0;
        int gap = 0;

        for (const auto &e : *this)
        {
            if (e == '=')
                score += params.match, gap = 0;
            else if (e == 'X')
                score += params.mismatch, gap = 0;
            else
            {
                if (gap == 0)
                    score -= params.gap_opening, gap++;
                score -= params.gap_extension;
            }
        }
        return score;
    }
};

using Cigars = std::vector<Cigar>;

struct nw_t
{
    int score{};
    Cigar cigar{};
    bool zdropped = false;

    nw_t() = default;
    explicit nw_t(const nw_t &nw) = default;
    explicit nw_t(nw_t &&nw) = default;
    nw_t &operator=(const nw_t &nw) = default;
    nw_t &operator=(nw_t &&nw) = default;
    ~nw_t() = default;

    nw_t(int s, const Cigar &c) : score(s), cigar(c) {}
};

/********** Set / Sequence **********/

using Sequence = std::string;
using Set = std::vector<Sequence>;
using Sets = std::vector<Set>;
using CompressedSequence = std::vector<uint8_t>;
using CompressedSequences = std::vector<uint8_t>;
using CompressedSet = std::vector<CompressedSequence>;

/**
 * @brief Returns number of unique pair that can be made from set of sequences.
 *
 * @param set set of sequence to pair
 * @return size_t
 */
inline size_t count_unique_pair(const Set &set)
{
    return (set.size() * (set.size() - 1)) / 2;
}

inline size_t count_unique_pair(const Sets &sets)
{
    size_t res = 0;
    for (const auto &set : sets)
        res += count_unique_pair(set);
    return res;
}

inline size_t count_compute_load(const Set &set)
{
    size_t compute_load = 0;
    for (size_t i = 0; i < set.size() - 1; i++)
        for (size_t j = i + 1; j < set.size(); j++)
            compute_load += set[i].size() + set[j].size() - 1;

    return compute_load;
}

inline size_t count_compute_load(const Sets &sets)
{
    size_t compute_load = 0;
    for (const auto &set : sets)
        compute_load += count_compute_load(set);

    return compute_load;
}

/********** nucleotide to ksw2 encoding **********/

/**
 * @brief Encode from ACTG/actg to 0123, undefined outside good values
 *
 * @param c nucleotide
 */
constexpr inline void encode_base(char &c)
{
    c >>= 1;
    constexpr uint8_t mask = 0b11;
    c &= mask;
}

/**
 * @brief Encode from ACTG to 0123 all element in collections.
 *
 * @param S
 */
inline void encode_base(auto &S)
{
    for (auto &s : S)
        encode_base(s);
}

inline Set encode_set(Set &&S)
{
    for (auto &s : S)
        encode_base(s);
    return S;
}

inline Sets encode_sets(Sets &&S)
{
    for (auto &s : S)
        encode_base(s);
    return S;
}

template <typename T, typename F>
constexpr inline auto operator|(T &&t, F f)
{
    return f(std::forward<T>(t));
}

void dump_to_file(const std::filesystem::path &filename, const auto &Container, const auto &Accessor)
{
    std::ofstream file(filename);

    printf("Writing %s\n", filename.c_str());

    for (const auto &e : Container)
        file << Accessor(e) << '\n';
}

[[noreturn]] static inline void exit(std::string message)
{
    printf("%s\n", message.c_str());
    std::exit(EXIT_FAILURE);
}

template <typename Container>
static constexpr inline auto resize(size_t i)
{
    return [i](Container &&s) -> Container
    {
        if (s.size() > i)
            s.resize(i);
        return s;
    };
}

static inline auto print_sets_size(const std::string &str)
{
    return [str](Sets &&sets)
    {
        printf(("  " + str + ": %lu\n").c_str(), sets.size());
        return sets;
    };
}

static inline auto print_set_size(const std::string &str)
{
    return [str](Set &&set)
    {
        printf(("  " + str + ": %lu\n").c_str(), set.size());
        return set;
    };
}

static inline auto sum_integers(auto i)
{
    static_assert(std::is_integral<decltype(i)>::value, "Integral type required !");
    return i * (i - 1) / 2;
}

static inline size_t triangular_index(size_t i, size_t j, size_t n)
{
    return sum_integers(n) - sum_integers(n - i) + j - i - 1;
}

#endif /* AD18B383_F97E_4512_98DA_46CE2947ACDD */
