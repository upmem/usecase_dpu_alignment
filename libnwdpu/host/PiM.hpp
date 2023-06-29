#ifndef CF5EAFF9_4B2B_4887_B1B9_EB3B3608047F
#define CF5EAFF9_4B2B_4887_B1B9_EB3B3608047F

#include <chrono>
#include <numeric>
#include <thread>

#include "Rank.hpp"

template <class Algo>
class PiM
{
    std::filesystem::path bin_path;
    std::vector<Rank<Algo>> m_ranks;

public:
    explicit PiM(const std::filesystem::path &filename, size_t n) : bin_path(filename), m_ranks(n)
    {
#pragma omp parallel for num_threads(4)
        for (auto &r : m_ranks)
            r.init(filename);
    }

    template <typename T>
    void send_all(T &data, const std::string &symbol)
    {
        for (auto &r : m_ranks)
            r.send_all(data, symbol);
    }

    Rank<Algo> &get_free_rank()
    {
        using namespace std::chrono_literals;

        while (true)
        {
            for (auto &r : m_ranks)
                if (r.is_available())
                    return r.alot();
            std::this_thread::sleep_for(1ms);
        }
    }

    void sync()
    {
        for (auto &r : m_ranks)
            dpu_sync(r.get());
    }

    void Print()
    {
        auto n_dpu = std::accumulate(m_ranks.begin(), m_ranks.end(), 0LU, [](size_t i, const auto &e)
                                     { return i + e.size(); });
        printf("PiM Accelerator with %lu ranks (%lu dpus).\n", m_ranks.size(), n_dpu);
    }
};

#endif /* CF5EAFF9_4B2B_4887_B1B9_EB3B3608047F */
