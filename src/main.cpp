/*
 * Copyright 2022 - UPMEM
 */

#include <fstream>
#include <yaml-cpp/yaml.h>

#include "../libnwdpu/host/dpu_common.hpp"
#include "fasta.hpp"

auto read_parameters(const std::filesystem::path &filename)
{
    auto config = YAML::LoadFile(filename.native());
    auto path = config["dataset"].as<std::string>();
    auto sets_number = config["sets_number"].as<uint32_t>();
    auto dpus_number = config["dpus_number"].as<uint32_t>();
    auto params = config["nw_params"];

    return std::tuple{
        path,
        sets_number,
        NW_Parameters{params["match"].as<int>(),
                      params["mismatch"].as<int>(),
                      params["gap_opening"].as<int>(),
                      params["gap_extension"].as<int>(), 128},
        dpus_number};
}

int main()
{
    const auto home = std::filesystem::canonical("/proc/self/exe").parent_path();

    auto [dataset_path, nsets, nw_parameters, ndpus] = read_parameters(home / "params.yaml");

    printf("DPU mode:\n"
           "  forcing width to 128.\n"
           "  using %u dpus.\n\n",
           ndpus);
    nw_parameters.width = 128;
    nw_parameters.print();

    printf("Dataset:\n");
    Timer load_time{};
    auto dataset = read_set_fasta(home / dataset_path) |
                   print_size("max") |
                   resize(nsets) |
                   print_size("use") |
                   encode_sets;
    load_time.print("  ");

    Timer compute_time{};
    auto alignments = dpu_pipeline("./libnwdpu/dpu/nw_affine", nw_parameters, ndpus, dataset);
    compute_time.print("  ");

    dump_to_file("scores.txt", alignments, [](const auto &e)
                 { return e.score; });

    dump_to_file("cigars.txt", alignments, [](const auto &e)
                 { return e.cigar; });

    return 0;
}
