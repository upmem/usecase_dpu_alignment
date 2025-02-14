/*
 * Copyright 2022 - UPMEM
 */

#include <fstream>

#include "../libnwdpu/host/dpu_common.hpp"
#include "timeline.hpp"
#include <vector>
#include <variant>
#include <functional>

#include "parameters.hpp"
template <class... Ts>
struct overloaded : Ts...
{
    using Ts::operator()...;
};
template <class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;
using AlignmentResult = std::variant<std::vector<NwType>, std::vector<int>>;

void write_to_file(AlignmentResult &alignments)
{
    std::visit(overloaded{[](const std::vector<NwType> &vec)
                          {
                              dump_to_file("scores.txt", vec, [](const auto &e)
                                           { return e.score; });
                              dump_to_file("cigars.txt", vec, [](const auto &e)
                                           { return e.cigar; });
                          },
                          [](const std::vector<int> &vec)
                          {
                              dump_to_file("scores.txt", vec, [](const auto &e)
                                           { return e; });
                          }},
               alignments);
}

int main(int argc, char **argv)
{
    auto [dataset_path, nsets, nw_parameters, ranks, app_mode] = populate_parameters(argc, argv);

    Timeline timeline{"log_times.csv"};

    printf("DPU ranks: %u\n\n", ranks);
    nw_parameters.Print();

    printf("Dataset:\n");
    Timer load_time{};
    auto dataset = read_set_fasta(dataset_path) |
                   print_size<Sets>("  max: ") |
                   resize<Sets>(nsets) |
                   print_size<Sets>("  use: ") |
                   encode<Sets>;
    load_time.Print("  ");

    timeline.mark("Initialization");

    AlignmentResult alignments;

    Timer compute_time{};
    switch (app_mode)
    {
    case AppMode::Set:
    {
        alignments = dpu_cigar_pipeline("./libnwdpu/dpu/nw_affine", nw_parameters, ranks, dataset);
        break;
    }
    case AppMode::Pair:
        break;
    case AppMode::All:
    {
        printf("All against all mode not implemented yet, use dpu_16S\n");
        break;
    }
    default:
        printf("Unknown application mode\n");
        break;
    }

    compute_time.Print("  ");

    timeline.mark("Alignement");

    if (std::holds_alternative<std::vector<NwType>>(alignments))
    {
        printf("Returned vector of NwType\n");
        printf("Scores and CIGARs were computed\n");
    }
    else if (std::holds_alternative<std::vector<int>>(alignments))
    {
        printf("Returned vector of int\n");
        printf("Only scores were computed\n");
    }

    write_to_file(alignments);

    return 0;
}
