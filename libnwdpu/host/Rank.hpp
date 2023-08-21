#ifndef BFFF1F0E_2A88_4908_A231_706B5C411C97
#define BFFF1F0E_2A88_4908_A231_706B5C411C97

#include <filesystem>

#include "../../src/types.hpp"

extern "C"
{
#include <dpu.h>
}

dpu_set_t init_rank(const std::filesystem::path &filename)
{
    dpu_set_t dpu{};
    if (!std::filesystem::exists(filename))
        exit("File " + filename.string() + " not found.");

    DPU_ASSERT(dpu_alloc_ranks(1, NULL, &dpu));
    DPU_ASSERT(dpu_load(dpu, filename.c_str(), NULL));
    return dpu;
}

size_t rank_size(dpu_set_t rank)
{
    uint32_t size = 0;
    dpu_get_nr_dpus(rank, &size);
    return size;
}

std::vector<size_t> rank_sizes(dpu_set_t ranks)
{
    dpu_set_t rank{};
    std::vector<size_t> sizes;

    DPU_RANK_FOREACH(ranks, rank)
    {
        sizes.push_back(rank_size(rank));
    }

    return sizes;
}

template <class App>
class Rank
{
    dpu_set_t m_rank{};
    size_t m_size{};
    bool m_valid = false;
    bool m_available = true;

public:
    App algo{};

    inline void init(const std::filesystem::path &filename)
    {
        m_rank = init_rank(filename);
        m_valid = true;
        m_size = rank_size(m_rank);
        algo.init(m_size);
    }

    size_t size() const { return m_size; }
    bool is_available() const { return m_available; }
    dpu_set_t &get() { return m_rank; }
    Rank &alot()
    {
        m_available = false;
        return *this;
    }
    void done() { m_available = true; }

    static dpu_error_t rank_done([[maybe_unused]] dpu_set_t _, [[maybe_unused]] uint32_t id, void *_arg)
    {
        auto &rank = *static_cast<Rank<App> *>(_arg);
        rank.done();
        return DPU_OK;
    };

    void launch() { DPU_ASSERT(dpu_launch(m_rank, DPU_ASYNCHRONOUS)); }
    void send() { algo.send(*this); }

    template <typename T>
    void send_all(T &data, const std::string &symbol)
    {
        if constexpr (requires(T t) { t.size(); t.data(); })
            DPU_ASSERT(dpu_broadcast_to(m_rank, symbol.c_str(), 0, data.data(), data.size() * sizeof(decltype(data[0])), DPU_XFER_ASYNC));
        else
            DPU_ASSERT(dpu_broadcast_to(m_rank, symbol.c_str(), 0, &data, sizeof(T), DPU_XFER_ASYNC));
    }

    void gather() { algo.gather(*this); }
    void post()
    {
        algo.post(*this);
        DPU_ASSERT(dpu_callback(m_rank, Rank<App>::rank_done, this, DPU_CALLBACK_ASYNC));
    };

    Rank() = default;
    Rank(const Rank &) = delete;
    Rank(Rank &&) = delete;
    Rank &operator=(const Rank &) = delete;
    Rank &operator=(Rank &&) = delete;
    ~Rank()
    {
        if (m_valid)
            dpu_free(m_rank);
    }
};

#endif /* BFFF1F0E_2A88_4908_A231_706B5C411C97 */
