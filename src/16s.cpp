/*
 * Copyright 2022 - UPMEM
 */

#include <cassert>
#include <fstream>
#include <yaml-cpp/yaml.h>

#include "../libnwdpu/host/dpu_common.hpp"
#include "fasta.hpp"

auto read_parameters(const std::filesystem::path &filename)
{
    auto config = YAML::LoadFile(filename.native());
    auto dataset = config["dataset"].as<std::string>();
    auto dpus_number = config["dpus_number"].as<uint32_t>();

    const auto home = std::filesystem::canonical("/proc/self/exe").parent_path();

    NW_Parameters nwp(
        config["nw_params"]["match"].as<int>(),
        config["nw_params"]["mismatch"].as<int>(),
        config["nw_params"]["gap_opening"].as<int>(),
        config["nw_params"]["gap_extension"].as<int>(),
        config["nw_params"]["width"].as<int>());

    return std::tuple{home / dataset, nwp, dpus_number};
}

int main()
{
    auto [dataset_path, params, ndpu] = read_parameters("./16s.yaml");

    printf("DPU mode:\n"
           "  forcing width to 128.\n"
           "  using %u dpus.\n\n",
           ndpu);
    params.width = 128;
    params.print();

    printf("Dataset:\n");
    auto dataset = read_seq_fasta(dataset_path) |
                   print_set_size("size") |
                   encode_set;

    Timer compute_time{};
    auto alignments = dpu_pipeline_16s("./libnwdpu/dpu/nw_16s", params, ndpu, dataset);
    compute_time.print("  ");

    dump_to_file("scores.txt", alignments, [](const auto &e)
                 { return e; });

    return 0;
}
