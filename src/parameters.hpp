#ifndef B31A7004_1AB6_4DC0_A536_DAB73DAC2F8E
#define B31A7004_1AB6_4DC0_A536_DAB73DAC2F8E

#include <filesystem>
#include <yaml-cpp/yaml.h>

#include "cxxopts.hpp"
#include "fasta.hpp"

enum class AppMode
{
    Set,
    Pair,
    All
};

inline std::ostream &operator<<(std::ostream &os, const AppMode &mode)
{
    switch (mode)
    {
    case AppMode::Set:
        os << "set";
        break;
    case AppMode::Pair:
        os << "pair";
        break;
    case AppMode::All:
        os << "all";
        break;
    }
    return os;
}

inline std::istream &operator>>(std::istream &is, AppMode &mode)
{
    std::string token;
    is >> token;
    if (token == "set")
        mode = AppMode::Set;
    else if (token == "pair")
        mode = AppMode::Pair;
    else if (token == "all")
        mode = AppMode::All;
    else
        throw std::runtime_error("Invalid app mode");
    return is;
}

inline std::string to_string(const AppMode &mode)
{
    switch (mode)
    {
    case AppMode::Set:
        return "set";
    case AppMode::Pair:
        return "pair";
    case AppMode::All:
        return "all";
    }
    throw std::runtime_error("Invalid app mode");
}

auto read_parameters(const std::filesystem::path &filename)
{
    auto config = YAML::LoadFile(filename.native());
    auto path = config["dataset"].as<std::string>();
    auto sets_number = config["sets_number"].as<uint32_t>();
    auto ranks = config["ranks"].as<uint32_t>();
    auto params = config["nw_params"];
    auto app_mode_str = config["app_mode"].as<std::string>();
    AppMode app_mode;

    if (app_mode_str == "set")
        app_mode = AppMode::Set;
    else if (app_mode_str == "pair")
        app_mode = AppMode::Pair;
    else if (app_mode_str == "all")
        app_mode = AppMode::All;
    else
        throw std::runtime_error("Invalid app mode in config file");

    return std::tuple{
        path,
        sets_number,
        NwParameters{params["match"].as<int32_t>(),
                     params["mismatch"].as<int32_t>(),
                     params["gap_opening"].as<int32_t>(),
                     params["gap_extension"].as<int32_t>(),
                     128},
        ranks,
        app_mode};
}

cxxopts::ParseResult parse_command_line(int argc, char **argv)
{
    cxxopts::Options options("sets", "Needleman-Wunsch algorithm on sets of sequences");
    options.add_options()(
        "c,config", "Path to the YAML configuration file (optional, superseded if other options are provided)", cxxopts::value<std::string>())(
        "d,dataset", "Path to the dataset file", cxxopts::value<std::string>())(
        "n,sets_number", "Number of sets to process", cxxopts::value<uint32_t>())(
        "r,ranks", "Number of ranks to use", cxxopts::value<uint32_t>())(
        "m,match", "Match score", cxxopts::value<int32_t>())(
        "x,mismatch", "Mismatch score", cxxopts::value<int32_t>())(
        "g,gap_opening", "Gap opening score", cxxopts::value<int32_t>())(
        "e,gap_extension", "Gap extension score", cxxopts::value<int32_t>())(
        "a,app_mode", "Application mode (set, pair, all)", cxxopts::value<AppMode>());

    options.add_options()("h,help", "Print usage");

    auto result = options.parse(argc, argv);

    if (result.count("help"))
    {
        printf("%s\n", options.help().c_str());
        exit(0);
    }

    if (!result.count("config") && (!result.count("dataset") || !result.count("sets_number") || !result.count("ranks") ||
                                    !result.count("match") || !result.count("mismatch") || !result.count("gap_opening") || !result.count("gap_extension") || !result.count("app_mode")))
    {
        printf("%s\n", options.help().c_str());
        exit(0);
    }

    return result;
}

template <typename T>
void update_parameter(const cxxopts::ParseResult &result, const std::string &option, T &parameter)
{
    if (result.count(option))
    {
        parameter = result[option].as<T>();
    }
}

auto populate_parameters(int argc, char **argv)
{
    auto result = parse_command_line(argc, argv);

    std::string path{};
    uint32_t sets_number{};
    uint32_t ranks{};
    NwParameters nw_parameters{0, 0, 0, 0, 128};
    AppMode app_mode{AppMode::Set};

    if (result.count("config") > 0)
    {
        std::tie(path, sets_number, nw_parameters, ranks, app_mode) = read_parameters(result["config"].as<std::string>());
    }

    update_parameter(result, "dataset", path);
    update_parameter(result, "sets_number", sets_number);
    update_parameter(result, "ranks", ranks);
    update_parameter(result, "match", nw_parameters.match);
    update_parameter(result, "mismatch", nw_parameters.mismatch);
    update_parameter(result, "gap_opening", nw_parameters.gap_opening);
    update_parameter(result, "gap_extension", nw_parameters.gap_extension);
    update_parameter(result, "app_mode", app_mode);

    return std::tuple{
        path,
        sets_number,
        nw_parameters,
        ranks,
        app_mode};
}

#endif /* B31A7004_1AB6_4DC0_A536_DAB73DAC2F8E */
