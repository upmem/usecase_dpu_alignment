/*
 * Copyright 2022 - UPMEM
 */

#include <cassert>
#include <yaml-cpp/yaml.h>

#include "../libnwdpu/host/dpu_common.hpp"
#include "fasta.hpp"
#include "../timeline.hpp"

auto read_parameters(const std::filesystem::path &filename)
{
    auto config = YAML::LoadFile(filename.native());
    auto dataset = config["dataset"].as<std::string>();
    auto ranks = config["ranks"].as<uint32_t>();
    auto params = config["nw_params"];

    const auto home = std::filesystem::canonical("/proc/self/exe").parent_path();

    return std::tuple{
        home / dataset,
        NwParameters{params["match"].as<int32_t>(),
                     params["mismatch"].as<int32_t>(),
                     params["gap_opening"].as<int32_t>(),
                     params["gap_extension"].as<int32_t>(),
                     128},
        ranks};
}

int main()
{
    auto [dataset_path, params, ranks] = read_parameters("./16s.yaml");
    Timeline timeline{"sets_time.csv"};

    printf("DPU mode:\n"
           "  forcing width to 128.\n"
           "  using %u ranks.\n\n",
           ranks);
    params.Print();

    printf("Dataset:\n");
    auto dataset = read_seq_fasta(dataset_path) |
                   print_size<Set>("  size: ") |
                   encode<Set>;

    timeline.mark("Initialization");
    Timer compute_time{};
    auto alignments = dpu_16s_pipeline("./libnwdpu/dpu/nw_16s", params, ranks, dataset);
    compute_time.Print("  ");
    timeline.mark("Alignement");

    dump_to_file("scores.txt", alignments, [](const auto &e)
                 { return e; });

    return 0;
}
