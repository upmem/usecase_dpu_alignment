#ifndef C9A0FB00_AA91_4E0A_B093_5EB1D2AA8D50
#define C9A0FB00_AA91_4E0A_B093_5EB1D2AA8D50

extern "C"
{
#include <time.h>
}

#include <fstream>
#include <string>

class Timeline
{
public:
    Timeline(const std::string &output_filename) : file(output_filename) { file << get_current_timestamp() << ",\n"; }

    void mark(const std::string &label) { file << get_current_timestamp() << ',' << label << '\n'; }

private:
    std::ofstream file;

    std::string get_current_timestamp()
    {
        struct timespec now;
        timespec_get(&now, TIME_UTC);
        return std::to_string(now.tv_sec) + '.' + std::to_string(now.tv_nsec);
    }
};

#endif /* C9A0FB00_AA91_4E0A_B093_5EB1D2AA8D50 */
