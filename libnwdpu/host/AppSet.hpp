#ifndef B264D0F9_04DC_4735_A7FD_3F0F386D53EE
#define B264D0F9_04DC_4735_A7FD_3F0F386D53EE

#include <algorithm>
#include <cassert>
#include <span>

#include "dpu_common.hpp"
#include "Rank.hpp"

struct SortedMap
{
    size_t index{};
    size_t load{};
    size_t dpu{};
    size_t offset{};
};

constexpr size_t max_dpu_load = 3000000;

inline auto sorted_map(const Sets &data)
{
    std::vector<SortedMap> index{data.size()};

    size_t offset = 0;
    for (size_t i = 0; i < index.size(); i++)
    {
        index[i] = {i, count_compute_load(data[i]), 0, offset};
        offset += sum_integers(data[i].size());
    }

    std::ranges::sort(index, [](const auto &a, const auto &b)
                      { return a.load > b.load; });

    return index;
}

size_t count_compute_load(const std::span<SortedMap> &index)
{
    size_t load = 0;
    for (const auto &e : index)
        load += e.load;
    return load;
}

auto dpu_to_cpu(std::vector<std::vector<NwType>> &res, const NwInputCigar &input, const NwCigarOutput &output, size_t mi)
{
    const auto &meta = input.metadata;
    res.resize(meta.number_of_sets);

    auto comb = [](auto n)
    { return n * (n - 1) / 2; };

    int dpu_offset = 0;

    for (size_t set_id = 0; set_id < meta.number_of_sets; set_id++)
    {
        int set_offset = 0;
        res[set_id].resize(comb(meta.set_sizes[set_id]));
        for (int i = 0; i < meta.set_sizes[set_id]; i++)
            for (int j = i + 1; j < meta.set_sizes[set_id]; j++)
            {
                res[set_id][set_offset].score = output.scores[dpu_offset];
                res[set_id][set_offset].dpu_offset = dpu_offset;
                res[set_id][set_offset].mi = mi;
                dpu_offset++;
                set_offset++;
            }
    }

    return res;
}

class AppSet
{
public:
    std::vector<NwInputCigar> inputs{};
    std::vector<NwCigarOutput> outputs{};
    std::vector<std::vector<char>> cigars{};
    std::span<SortedMap> index{};
    std::span<NwType> result{};
    size_t cigar_size{};

    inline void init(size_t size)
    {
        inputs.resize(size);
        outputs.resize(size);
        cigars.resize(size);
    }

    void send(Rank<AppSet> &rank)
    {
        dpu_set_t dpu{};
        uint32_t each_dpu = 0;

        DPU_FOREACH(rank.get(), dpu, each_dpu)
        {
            DPU_ASSERT(dpu_prepare_xfer(dpu, inputs[each_dpu].sequences.data()));
        }
        DPU_ASSERT(dpu_push_xfer(rank.get(), DPU_XFER_TO_DPU, "sequences", 0,
                                 SCORE_MAX_SEQUENCES_TOTAL_SIZE, DPU_XFER_ASYNC));

        DPU_FOREACH(rank.get(), dpu, each_dpu)
        {
            DPU_ASSERT(dpu_prepare_xfer(dpu, &inputs[each_dpu].metadata));
        }
        DPU_ASSERT(dpu_push_xfer(rank.get(), DPU_XFER_TO_DPU, "metadata", 0,
                                 sizeof(NwMetadataDPU), DPU_XFER_ASYNC));

        DPU_FOREACH(rank.get(), dpu, each_dpu)
        {
            DPU_ASSERT(dpu_prepare_xfer(dpu, inputs[each_dpu].cigar_indexes.data()));
        }
        DPU_ASSERT(dpu_push_xfer(rank.get(), DPU_XFER_TO_DPU, "cigar_indexes", 0,
                                 METADATA_MAX_NUMBER_OF_SCORES * sizeof(uint32_t), DPU_XFER_ASYNC));
    }

    void gather(Rank<AppSet> &rank)
    {
        dpu_set_t dpu{};
        uint32_t each_dpu = 0;

        for (auto &cigar : cigars)
            cigar.resize(cigar_size);

        DPU_FOREACH(rank.get(), dpu, each_dpu)
        {
            DPU_ASSERT(dpu_prepare_xfer(dpu, &outputs[each_dpu]));
        }
        DPU_ASSERT(dpu_push_xfer(rank.get(), DPU_XFER_FROM_DPU, "output", 0,
                                 sizeof(NwCigarOutput), DPU_XFER_ASYNC));
        DPU_FOREACH(rank.get(), dpu, each_dpu)
        {
            DPU_ASSERT(dpu_prepare_xfer(dpu, cigars[each_dpu].data()));
        }
        DPU_ASSERT(dpu_push_xfer(rank.get(), DPU_XFER_FROM_DPU, "cigars", 0,
                                 cigar_size, DPU_XFER_ASYNC));
    }

    void post(Rank<AppSet> &rank)
    {
        DPU_ASSERT(dpu_callback(rank.get(), rank_postprocess, this, DPU_CALLBACK_ASYNC));
    }

    static dpu_error_t rank_postprocess([[maybe_unused]] dpu_set_t _, [[maybe_unused]] uint32_t id, void *_arg)
    {
        auto &rank = *static_cast<AppSet *>(_arg);
        const auto &inputs = rank.inputs;
        const auto &outputs = rank.outputs;
        auto &cpu_output = rank.result;
        const auto &cigars = rank.cigars;
        const auto &index = rank.index;

        auto size = inputs.size();

        std::array<std::vector<std::vector<NwType>>, 64> dpu_res{};
        // #pragma omp parallel for
        for (size_t i = 0; i < size; i++)
            dpu_to_cpu(dpu_res[i], inputs[i], outputs[i], i);

        for (const auto &[i, load, d, off] : index)
        {
            auto &set_res = dpu_res[d].front();
            for (size_t r = 0; r < set_res.size(); r++)
            {
                cpu_output[off + r] = set_res[r];
                cpu_output[off + r].cigar.resize(outputs[set_res[r].mi].lengths[set_res[r].dpu_offset]);
                cpu_output[off + r].cigar.assign(&cigars[set_res[r].mi][inputs[set_res[r].mi].cigar_indexes[set_res[r].dpu_offset]], outputs[set_res[r].mi].lengths[set_res[r].dpu_offset]);
            }
            dpu_res[d].erase(dpu_res[d].begin());
        }

        return DPU_OK;
    }

    static auto take_load(std::span<SortedMap> &sp, size_t load, size_t n_dpu)
    {
        if (n_dpu == 0)
            exit("Rank size is 0 !\n");

        int i = 0;
        size_t tot_load = 0;

        for (const auto &sm : sp)
        {
            if (tot_load >= load && i >= static_cast<int>(n_dpu))
                break;

            tot_load += sm.load;
            i++;
        }

        if (i > static_cast<int>(n_dpu))
            i = i - static_cast<int>(i % n_dpu);

        auto b = sp.begin();

        sp = sp.last(sp.size() - i);

        return std::span<SortedMap>(b, b + i);
    }

    void get_bucket(std::span<SortedMap> &index_span)
    {
        static size_t thres = 16000000;

        auto n_dpu = inputs.size();

        auto threshold = thres;
        thres = thres - (thres / 45);
        auto rl = threshold * n_dpu;

        index = take_load(index_span, rl, n_dpu);
    }

    static inline auto cpu_to_dpu(const Sets &sets, NwInputCigar &dpu_input)
    {
        assert(sets.size() <= SCORE_METADATA_MAX_NUMBER_OF_SET &&
               "Too many sets for DPU!\n");

        auto &meta = dpu_input.metadata;
        meta.number_of_sets = static_cast<uint32_t>(sets.size());

        uint32_t idx = 0;
        uint32_t seq_idx = 0;
        uint32_t cigar_index = 0;
        uint32_t cigar_offset = 0;

        for (size_t set_id = 0; set_id < sets.size(); set_id++)
            meta.set_sizes[set_id] = static_cast<uint8_t>(sets[set_id].size());

        for (const auto &set : sets)
        {
            for (const auto &seq : set)
            {
                auto cseq = compress_sequence(seq);
                meta.lengths[seq_idx] = static_cast<uint16_t>(seq.size());
                meta.indexes[seq_idx++] = idx;
                idx += cseq.size();
                push_back(dpu_input.sequences, cseq);
            }
            assert(dpu_input.sequences.size() < SCORE_MAX_SEQUENCES_TOTAL_SIZE &&
                   "dpu sequence buffer overflow!\n");

            for (size_t i = 0; i < set.size() - 1; i++)
                for (size_t j = i + 1; j < set.size(); j++)
                {
                    dpu_input.cigar_indexes[cigar_offset] = cigar_index;
                    cigar_offset++;
                    auto max_cigar_size = set[i].size() + set[j].size();

                    assert(set[i].size() + set[j].size() < UINT16_MAX && "cigar is to big for uint16_t\n");

                    max_cigar_size += (8 - (max_cigar_size % 8));
                    cigar_index += max_cigar_size;

                    assert(cigar_index < MAX_CIGAR_SIZE && "not enough space for cigar on dpu!\n");
                    assert(cigar_offset <= METADATA_MAX_NUMBER_OF_SCORES && "too much cigar to compute!\n");
                }
        }

        return size_t{cigar_index};
    }

    static auto bucket_sets(const Sets &data, auto &index, size_t n)
    {
        std::vector<Sets> dpu_sets(n);
        std::vector<size_t> dpu_loads(n);

        for (auto &[i, load, d, off] : index)
        {
            auto min_index = std::distance(dpu_loads.begin(), std::ranges::min_element(dpu_loads));
            dpu_sets[min_index].push_back(data[i]);
            dpu_loads[min_index] += load;
            d = min_index;
        }

        size_t s = 0;
        for (const auto &set : dpu_sets)
            s += set.size();

        return dpu_sets;
    }

    void to_dpu_format(const Sets &data, const NwParameters &p)
    {
        auto dpu_sets = bucket_sets(data, index, inputs.size());

        cigar_size = 0;

        for (size_t i = 0; i < inputs.size(); i++)
        {
            inputs[i].metadata.match = p.match;
            inputs[i].metadata.mismatch = p.mismatch;
            inputs[i].metadata.gap_opening = p.gap_opening;
            inputs[i].metadata.gap_extension = p.gap_extension;
            inputs[i].sequences.resize(0);
            inputs[i].sequences.reserve(SCORE_MAX_SEQUENCES_TOTAL_SIZE);
            inputs[i].cigar_indexes.resize(METADATA_MAX_NUMBER_OF_SCORES);

            cigar_size = std::max(cpu_to_dpu(dpu_sets[i], inputs[i]), cigar_size);
        }
    }
};

#endif /* B264D0F9_04DC_4735_A7FD_3F0F386D53EE */
