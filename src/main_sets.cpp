/*
 * Copyright 2022 - UPMEM
 */

#include <fstream>
#include <yaml-cpp/yaml.h>

#include "../libnwdpu/host/dpu_common.hpp"
#include "fasta.hpp"
#include "../timeline.hpp"

auto read_parameters(const std::filesystem::path &filename)
{
    auto config = YAML::LoadFile(filename.native());
    auto path = config["dataset"].as<std::string>();
    auto sets_number = config["sets_number"].as<uint32_t>();
    auto ranks = config["ranks"].as<uint32_t>();
    auto params = config["nw_params"];

    return std::tuple{
        path,
        sets_number,
        NwParameters{params["match"].as<int32_t>(),
                     params["mismatch"].as<int32_t>(),
                     params["gap_opening"].as<int32_t>(),
                     params["gap_extension"].as<int32_t>(),
                     128},
        ranks};
}

int main()
{
    const auto home = std::filesystem::canonical("/proc/self/exe").parent_path();

    auto [dataset_path, nsets, nw_parameters, ranks] = read_parameters(home / "sets.yaml");

    Timeline timeline{"sets_time.csv"};

    printf("DPU mode:\n"
           "  Forcing width to 128.\n"
           "  Asking for %u ranks.\n\n",
           ranks);
    nw_parameters.Print();

    printf("Dataset:\n");
    Timer load_time{};
    auto dataset = read_set_fasta(home / dataset_path) |
                   print_size<Sets>("  max: ") |
                   resize<Sets>(nsets) |
                   print_size<Sets>("  use: ") |
                   encode<Sets>;
    load_time.Print("  ");

    timeline.mark("Initialization");

    Timer compute_time{};
    auto alignments = dpu_cigar_pipeline("./libnwdpu/dpu/nw_affine", nw_parameters, ranks, dataset);
    compute_time.Print("  ");

    timeline.mark("Alignement");

    dump_to_file("scores.txt", alignments, [](const auto &e)
                 { return e.score; });

    /*dump_to_file("cigars.txt", alignments, [](const auto &e)
                 { return e.cigar; });*/

    return 0;
}
