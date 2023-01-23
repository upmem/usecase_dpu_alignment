/*
 * Copyright 2022 - UPMEM
 */

#ifndef AD18B383_F97E_4512_98DA_46CE2947ACDD
#define AD18B383_F97E_4512_98DA_46CE2947ACDD

#include <array>
#include <vector>

#include "../cdefs.h"
#include "timer.hpp"

static constexpr size_t NB_NUCLEOTIDE = 5;

using SimilarityMatrix = std::array<int8_t, NB_NUCLEOTIDE * NB_NUCLEOTIDE>;
constexpr SimilarityMatrix make_similarity_matrix(int pM, int pX)
{
    int8_t M = static_cast<int8_t>(pM);
    int8_t X = static_cast<int8_t>(pX);
    return {
        M, X, X, X, 0,
        X, M, X, X, 0,
        X, X, M, X, 0,
        X, X, X, M, 0,
        0, 0, 0, 0, 0};
}

struct NW_Parameters
{
    int match;
    int mismatch;
    int32_t gap_opening1;
    int32_t gap_extension1;
    int32_t gap_opening2;
    int32_t gap_extension2;
    int width;
    SimilarityMatrix smat;

    NW_Parameters() = default;
    NW_Parameters(const NW_Parameters &) = default;
    NW_Parameters(NW_Parameters &&) = default;
    NW_Parameters &operator=(const NW_Parameters &) = default;
    NW_Parameters &operator=(NW_Parameters &&) = default;
    ~NW_Parameters() = default;
    void print() const
    {
        printf("Alignment parameters:\n"
               "  match:         %d\n"
               "  mismatch:      %d\n"
               "  gap opening:   %d\n"
               "  gap extension: %d\n"
               "  gap opening:   %d\n"
               "  gap extension: %d\n"
               "  width:         %d\n\n",
               match, mismatch, gap_opening1, gap_extension1, gap_opening2, gap_extension2, width);
    }

    constexpr NW_Parameters(
        int _match,
        int _miss,
        int32_t _gapo1,
        int32_t _gape1,
        int32_t _gapo2,
        int32_t _gape2,
        int _width) : match(_match),
                      mismatch(_miss),
                      gap_opening1(_gapo1),
                      gap_extension1(_gape1),
                      gap_opening2(_gapo2),
                      gap_extension2(_gape2),
                      width(_width),
                      smat(make_similarity_matrix(_match, _miss))
    {
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
                    score -= params.gap_opening1, gap++;
                score -= params.gap_extension1;
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
using Sequences = std::vector<Sequence>;
using Set = Sequences;
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

inline Set encode_set(Set &S)
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

[[noreturn]] static inline void exit(std::string message)
{
    printf("%s\n", message.c_str());
    std::exit(EXIT_FAILURE);
}

#endif /* AD18B383_F97E_4512_98DA_46CE2947ACDD */
