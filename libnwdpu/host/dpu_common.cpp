/*
 * Copyright 2022 - UPMEM
 */

#include "dpu_common.hpp"
#include "PiM.hpp"
#include "AppSet.hpp"
#include "App16S.hpp"

extern "C"
{
#include <dpu.h>
}

std::vector<NwType> dpu_cigar_pipeline(std::filesystem::path dpu_bin_path, const NwParameters &p, size_t n_ranks, const Sets &sets)
{

    PiM<AppSet> accelerator(dpu_bin_path, n_ranks);
    accelerator.Print();

    auto index = sorted_map(sets);

    auto index_span = std::span<SortedMap>(index);
    std::vector<NwType> cpu_output(count_unique_pair(sets));

    size_t total_set = index_span.size();

    while (!index_span.empty())
    {
        auto &rank = accelerator.get_free_rank();
        rank.algo.get_bucket(index_span, total_set, n_ranks);
        rank.algo.result = cpu_output;
        rank.algo.to_dpu_format(sets, p);

        rank.send();
        rank.launch();
        rank.gather();
        rank.post();
    }

    accelerator.sync();

    // Add to dump DPU counters and analyse their individual workload.
    // dump_to_file("counters.txt", dpu_outputs, [](const auto &e) { return e.perf_counter; });

    return cpu_output;
}

auto Set_to_dpuSet(const Set &data, const NwParameters &params)
{
    NwInputScore dpu_input;

    dpu_input.sequences.reserve(SCORE_MAX_SEQUENCES_TOTAL_SIZE);
    dpu_input.metadata.match = params.match;
    dpu_input.metadata.mismatch = params.mismatch;
    dpu_input.metadata.gap_extension = params.gap_extension;
    dpu_input.metadata.gap_opening = params.gap_opening;

    size_t dpu_index = 0;
    size_t seq_id = 0;

    auto cset = compress_set(data);

    for (size_t i = 0; i < cset.size(); i++)
    {
        assert(seq_id < DPU_MAX_NUMBER_OF_SEQUENCES_MRAM && "Set is too big!\n");
        dpu_input.sequence_metadata.lengths[seq_id] = static_cast<uint16_t>(data[i].size());
        dpu_input.sequence_metadata.indexes[seq_id] = static_cast<uint32_t>(dpu_index);
        dpu_index += cset[i].size();
        seq_id++;
    }

    for (const auto &seq : cset)
        push_back(dpu_input.sequences, seq);

    assert(dpu_input.sequences.size() < SCORE_MAX_SEQUENCES_TOTAL_SIZE &&
           "dpu sequence buffer overflow!\n");

    return dpu_input;
}

std::vector<int> dpu_16s_pipeline(std::filesystem::path dpu_bin_path, const NwParameters &p, size_t n_ranks, const Set &set)
{

    PiM<App16S> accelerator(dpu_bin_path, n_ranks);
    accelerator.Print();

    auto dpu_dataset = Set_to_dpuSet(set, p);
    accelerator.send_all(dpu_dataset.sequences, "sequences");
    accelerator.send_all(dpu_dataset.metadata, "metadata");
    accelerator.send_all(dpu_dataset.sequence_metadata, "sequence_metadata");

    std::vector<int> cpu_output(sum_integers(set.size()));

    auto total_size = sum_integers(set.size());

    ComparisonMetadata meta{
        0,
        1,
        0,
        static_cast<uint32_t>(set.size())};

    while (total_size > 0)
    {
        auto thr = std::min(total_size / 80, 100000UL);
        auto i = std::max(thr, 512UL);
        i = std::min(i, total_size);
        total_size -= i;
        auto &rank = accelerator.get_free_rank();
        rank.algo.p_results = &cpu_output;
        rank.algo.get_bucket(meta, i);

        rank.send();
        rank.launch();
        rank.gather();
        rank.post();
    }

    accelerator.sync();

    return cpu_output;
}