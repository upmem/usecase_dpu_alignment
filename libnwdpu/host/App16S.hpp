#ifndef CACF6DF8_0DCE_4D9D_8EC6_5061ECA11076
#define CACF6DF8_0DCE_4D9D_8EC6_5061ECA11076

#include "dpu_common.hpp"
#include "Rank.hpp"

class App16S
{

public:
    std::vector<ComparisonMetadata> meta{};
    std::vector<NwScoreOutput> outputs{};
    std::vector<int> *p_results;

    inline void init(size_t size)
    {
        meta.resize(size);
        outputs.resize(size);
    }

    void send(Rank<App16S> &rank)
    {
        dpu_set_t dpu{};
        uint32_t each_dpu = 0;

        DPU_FOREACH(rank.get(), dpu, each_dpu)
        {
            DPU_ASSERT(dpu_prepare_xfer(dpu, &(meta[each_dpu])));
        }
        DPU_ASSERT(dpu_push_xfer(rank.get(), DPU_XFER_TO_DPU, "meta_index", 0,
                                 sizeof(ComparisonMetadata), DPU_XFER_ASYNC));
    }

    void gather(Rank<App16S> &rank)
    {
        dpu_set_t dpu{};
        uint32_t each_dpu = 0;

        auto byte_size = static_cast<size_t>(meta[0].count * 4 + 16);
        byte_size &= ~7;
        byte_size = std::min(byte_size, sizeof(NwScoreOutput));

        DPU_FOREACH(rank.get(), dpu, each_dpu)
        {
            DPU_ASSERT(dpu_prepare_xfer(dpu, &outputs[each_dpu]));
        }
        DPU_ASSERT(dpu_push_xfer(rank.get(), DPU_XFER_FROM_DPU, "output", 0, byte_size, DPU_XFER_ASYNC));
    }

    void post(Rank<App16S> &rank)
    {
        DPU_ASSERT(dpu_callback(rank.get(), rank_postprocess, this, DPU_CALLBACK_ASYNC));
    }

    static dpu_error_t rank_postprocess([[maybe_unused]] dpu_set_t _, [[maybe_unused]] uint32_t id, void *_arg)
    {
        auto &algo = *static_cast<App16S *>(_arg);

        for (size_t i = 0; i < algo.meta.size(); i++)
        {
            auto idx = triangular_index(algo.meta[i].start_row, algo.meta[i].start_col, algo.meta[i].size);

            for (size_t j = 0; j < algo.meta[i].count; j++)
            {
                algo.p_results->operator[](idx) = algo.outputs[i].scores[j];
                idx++;
            }
        }

        return DPU_OK;
    }

    static void update_meta(ComparisonMetadata &meta, int &rest)
    {
        for (uint32_t i = 0; i < meta.count; i++)
        {
            meta.start_col++;
            if (meta.start_col >= meta.size)
            {
                meta.start_row++;
                meta.start_col = meta.start_row + 1;
            }
        }

        if (rest-- == 0)
            meta.count--;
    }

    void get_bucket(ComparisonMetadata &new_meta, size_t i)
    {
        const auto nr_dpu = meta.size();
        auto mean = i / nr_dpu;
        auto rest = static_cast<int>(i % nr_dpu);

        new_meta.count = static_cast<uint32_t>(mean + (rest != 0 ? 1 : 0));
        rest--;

        for (auto &dpu : meta)
        {
            dpu = new_meta;
            update_meta(new_meta, rest);
        }
    }
};

#endif /* CACF6DF8_0DCE_4D9D_8EC6_5061ECA11076 */
