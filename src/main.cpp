/*
 * Copyright 2022 - UPMEM
 */

#include <yaml-cpp/yaml.h>

#include "../libnwdpu/host/dpu_common.hpp"
#include "fasta.hpp"

auto read_parameters(const std::filesystem::path &filename)
{
    auto config = YAML::LoadFile(filename.native());
    auto path = config["dataset"].as<std::string>();
    auto sets_number = config["sets_number"].as<uint32_t>();
    auto dpus_number = config["dpus_number"].as<uint32_t>();

    NW_Parameters nwp(
        config["nw_params"]["match"].as<int>(),
        config["nw_params"]["mismatch"].as<int>(),
        config["nw_params"]["gap_opening1"].as<int>(),
        config["nw_params"]["gap_extension1"].as<int>(),
        0, 0, 128);

    return std::tuple{path, sets_number, nwp, dpus_number};
}

static constexpr inline auto resize(size_t i)
{
    return [i](Sets &&s) -> Sets
    {
        if (s.size() > i)
            s.resize(i);
        return s;
    };
}

static inline auto print_dataset(Sets &&sets)
{
    printf("Dataset:\n"
           "   sets number: %lu\n",
           sets.size());
    return sets;
}

void print_cigar(NW_Parameters &params, const std::vector<nw_t> &cigar)
{
    printf("score: %d\n"
           "cigar: %.1000s\n"
           "cigar score: %d\n"
           "size:  %lu\n",
           cigar[1].score, cigar[1].cigar.c_str(), cigar[1].cigar.count_score(params), cigar[1].cigar.size());
}

int main()
{
    const auto home = std::filesystem::canonical("/proc/self/exe").parent_path();

    auto [dataset_path, nsets, nw_parameters, ndpus] = read_parameters(home / "params.yaml");

    printf("DPU mode:\n"
           "  forcing width to 128.\n"
           "  using %u dpus.\n",
           ndpus);
    nw_parameters.width = 128;
    nw_parameters.print();

    auto dataset = read_pacbio(home / dataset_path) |
                   resize(nsets) |
                   print_dataset |
                   encode_sets;

    Timer compute_time{};
    auto alignments = dpu_pipeline("./libnwdpu/dpu/nw_affine", nw_parameters, ndpus, dataset);
    compute_time.print();
    print_cigar(nw_parameters, alignments);

    return 0;
}
