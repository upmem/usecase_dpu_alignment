/*
 * Copyright 2022 - UPMEM
 */

#ifndef AD18B383_F97E_4512_98DA_46CE2947ACDD
#define AD18B383_F97E_4512_98DA_46CE2947ACDD

#include <array>
#include <concepts>
#include <filesystem>
#include <fstream>
#include <vector>

#include <immintrin.h>

#include "../cdefs.h"
#include "timer.hpp"

/**
 * @brief Needleman & Wunsch parameters
 *
 */
struct NwParameters
{
    /// @brief Public parameters
    int32_t match;         /// match bonus
    int32_t mismatch;      /// mismatch penalty
    int32_t gap_opening;   /// gap opening penalty
    int32_t gap_extension; /// gap extension penalty
    int32_t width;         /// band width

    /**
     * @brief Print Needleman & Wunsch parameters
     *
     */
    void Print() const
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

/**
 * @brief CIGAR Representation
 *
 */
struct Cigar : public std::string
{
    /**
     * @brief Computes score of CIGAR
     *
     * @param params Bonus/Penalties to use for score computation
     * @return Score
     */
    int CountScore(const NwParameters &params) const
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

/// @brief Type for a collection of CIGARs
using Cigars = std::vector<Cigar>;

/**
 * @brief Type for Needleman & Wunsch return values
 *
 */
struct NwType
{
    /// @brief Aggregate data
    int score{};   /// Score from alignment
    Cigar cigar{}; /// CIGAR of the alignment
    size_t dpu_offset{};
    size_t mi{};
};

/********** Set / Sequence **********/

using Sequence = std::string;                          /// Sequence type
using Set = std::vector<Sequence>;                     /// Set type
using Sets = std::vector<Set>;                         /// Sets type
using CompressedSequence = std::vector<uint8_t>;       /// Compressed sequence type
using CompressedSequences = std::vector<uint8_t>;      /// Compressed sequences type
using CompressedSet = std::vector<CompressedSequence>; /// Compressed set type

/**
 * @brief Sum of the first i numbers, starting at 0.
 *
 * @param i
 * @return auto
 */
static constexpr auto sum_integers(std::integral auto i)
{
    return i * (i - 1) / 2;
}

/**
 * @brief Round number to the next multiple of 8
 *
 */
static constexpr auto round_up8(std::integral auto n)
{
    return (n + 7) & ~7;
}

/**
 * @brief Returns the total number of unique pairs in a collection of set
 *
 * @param sets
 * @return size_t
 */
inline size_t count_unique_pair(const Sets &sets)
{
    size_t res = 0;
    for (const auto &set : sets)
        res += sum_integers(set.size());
    return res;
}

/**
 * @brief Estimate the number of cells to compute for given Set (assuming banded N&W)
 *
 * @param set
 * @return Load estimation
 */
inline size_t count_compute_load(const Set &set)
{
    size_t compute_load = 0;
    for (size_t i = 0; i < set.size() - 1; i++)
        for (size_t j = i + 1; j < set.size(); j++)
            compute_load += set[i].size() + set[j].size() - 1;

    return compute_load;
}

/**
 * @brief Estimate the total compute load of a collection of Set.
 *
 * @param sets
 * @return size_t
 */
inline size_t count_compute_load(const Sets &sets)
{
    size_t compute_load = 0;
    for (const auto &set : sets)
        compute_load += count_compute_load(set);

    return compute_load;
}

/********** nucleotide to ksw2 encoding **********/

/**
 * @brief Encode from ACTG/actg to 0123, undefined for other values
 *
 * @param c nucleotide
 */
template <typename C>
constexpr C &encode(C &&c)
{
    if constexpr (std::is_same_v<char, std::remove_reference_t<C>>)
    {
        c >>= 1;
        constexpr uint8_t mask = 0b11;
        c &= mask;

        return c;
    }
    else
    {
        for (auto &e : c)
            encode(e);

        return c;
    }
}

/**
 * @brief A small helper operator to chain operations
 *
 * @tparam T
 * @tparam F
 * @param t Data to apply the function on
 * @param f Next function to apply to t
 * @return constexpr auto
 */
template <typename T, typename F, typename std::enable_if_t<std::is_invocable<F, T>::value, bool> = true>
constexpr auto operator|(T &&t, F f)
{
    return f(std::forward<decltype(t)>(t));
}

/**
 * @brief Dump all the element from a container to a file. Accessor applied on element.
 *
 * @param filename
 * @param Container
 * @param Accessor
 */
void dump_to_file(const std::filesystem::path &filename, const auto &Container, const auto &Accessor)
{
    std::ofstream file(filename);

    printf("Writing %s\n", filename.c_str());

    for (const auto &e : Container)
        file << Accessor(e) << '\n';
}

/**
 * @brief Exit the program after printing a given error message
 *
 * @param message
 */
[[noreturn]] static inline void exit(std::string message)
{
    printf("%s\n", message.c_str());
    std::exit(EXIT_FAILURE);
}

/**
 * @brief Resize a container to a smaller size, nothing is done if size greater than current size.
 *
 * @tparam Container
 * @param i
 * @return constexpr auto
 */
template <typename Container>
static constexpr auto resize(size_t i)
{
    return [i](Container &&s) -> Container
    {
        if (s.size() > i)
            s.resize(i);
        return s;
    };
}

/**
 * @brief Return a ffunction printing the size of container with given string prefixed
 *
 * @tparam C
 * @param str
 * @return auto
 */
template <typename C>
static inline auto print_size(const std::string &str)
{
    return [str](C &&container)
    {
        printf((str + "%lu\n").c_str(), container.size());
        return container;
    };
}

/**
 * @brief Gives the index in a linear buffer of an upper triangular matrix index
 *
 * @param i i th row
 * @param j j th column
 * @param n matrix size
 * @return size_t
 */
static inline size_t triangular_index(size_t i, size_t j, size_t n)
{
    return sum_integers(n) - sum_integers(n - i) + j - i - 1;
}

/**
 * @brief Return the size a dpu buffer needs to contains the compressed representation of a sequence
 *
 * @param size
 * @return constexpr uint32_t
 */
constexpr uint32_t compressed_size(size_t size)
{
    auto compressed_size = static_cast<uint32_t>((size + 3) / 4);
    return round_up8(compressed_size);
}

/**
 * @brief Return the size a dpu buffer needs to contains the compressed representation of a set
 *
 * @param size
 * @return constexpr uint32_t
 */
constexpr uint32_t compressed_size(Set set)
{
    uint32_t total = 0;
    for (const auto &e : set)
        total += compressed_size(e.size());
    return total;
}

/**
 * @brief Return the size a dpu buffer needs to contains the compressed representation of a sets
 *
 * @param size
 * @return constexpr uint32_t
 */
constexpr uint32_t compressed_size(Sets sets)
{
    uint32_t total = 0;
    for (const auto &e : sets)
        total += compressed_size(e);
    return total;
}

/**
 * @brief Return the compressed representation a sequence
 *
 * @param seq
 * @return CompressedSequence
 */
inline CompressedSequence compress_sequence(const Sequence &seq)
{
    uint32_t csize = compressed_size(seq.size());
    CompressedSequence cseq(csize);

    for (size_t i = 0; i < csize; i += 2)
    {
        size_t seq_id = i * 4;

        auto c4n = _pext_u64(*(uint64_t *)(seq.data() + seq_id), 0x0303030303030303);
        *(uint16_t *)(cseq.data() + i) = c4n;
    }

    return cseq;
}

/**
 * @brief Returns a Set of the compressed sequences
 *
 * @param set
 * @return CompressedSet
 */
inline CompressedSet compress_set(const Set &set)
{
    CompressedSet cset(set.size());

#pragma omp parallel for
    for (size_t i = 0; i < set.size(); i++)
        cset[i] = compress_sequence(set[i]);

    return cset;
}

/**
 * @brief Concat a compressed sequence to an existing buffer
 *
 * @param cseqs
 * @param cseq
 */
inline void push_back(CompressedSequences &cseqs, const CompressedSequence &cseq)
{
    size_t begin = cseqs.size();
    cseqs.resize(begin + cseq.size());

    for (size_t i = 0; i < cseq.size(); i++)
        cseqs[begin + i] = cseq[i];
}

/**
 * @brief Concat a compressed sequence to an existing buffer
 *
 * @param cseqs
 * @param cseq
 */
inline uint32_t compressed_emplace(CompressedSequences &cseqs, const Sequence &seq)
{
    const uint32_t csize = compressed_size(seq.size());

    const size_t begin = cseqs.size();
    cseqs.resize(begin + csize);

    for (size_t i = 0; i < csize; i += 2)
    {
        size_t seq_id = i * 4;

        auto c4n = _pext_u64(*(uint64_t *)(seq.data() + seq_id), 0x0303030303030303);
        *(uint16_t *)(cseqs.data() + begin + i) = c4n;
    }

    return csize;
}

#endif /* AD18B383_F97E_4512_98DA_46CE2947ACDD */
