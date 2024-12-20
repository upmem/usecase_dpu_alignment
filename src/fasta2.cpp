#include <algorithm>
#include <fstream>
#include <sys/mman.h> // for mmap
#include <fcntl.h>    // for open
#include <unistd.h>   // for close

#include "fasta.hpp"

class MappedFile
{
    char *data;
    size_t size;

public:
    MappedFile() = delete;

    MappedFile(const std::string &filename)
    {
        int fd = open(filename.c_str(), O_RDONLY);
        if (fd == -1)
        {
            perror("Error opening file");
            exit(EXIT_FAILURE);
        }

        off_t s = lseek(fd, 0, SEEK_END);
        if (s == -1)
        {
            perror("Error getting file size");
            close(fd);
            exit(EXIT_FAILURE);
        }

        size = s;

        data = static_cast<char *>(mmap(nullptr, size, PROT_READ, MAP_PRIVATE, fd, 0));
        if (data == MAP_FAILED)
        {
            perror("Error mapping file");
            close(fd);
            exit(EXIT_FAILURE);
        }

        close(fd);
    }

    MappedFile(const MappedFile &) = delete;
    MappedFile(MappedFile &&) = delete;
    MappedFile &operator=(const MappedFile &) = delete;
    MappedFile &operator=(MappedFile &&) = delete;

    ~MappedFile()
    {
        if (munmap(data, size) == -1)
        {
            perror("Error unmapping file");
            exit(EXIT_FAILURE);
        }
    }

    inline const char *getData() const { return data; }
    inline size_t getSize() const { return size; }
};

std::vector<std::string> split_lines(const char *mapped_file, size_t file_size)
{
    std::vector<std::string> lines;
    const char *start = mapped_file;
    const char *end = mapped_file + file_size;
    while (start < end)
    {
        const char *line_end = std::find(start, end, '\n');
        lines.emplace_back(start, line_end);
        start = line_end + 1;
    }
    return lines;
}

static inline std::vector<std::string> check_lines_number(const std::vector<std::string> &lines)
{
    if ((lines.size() & 1) == 1)
        exit("Lines are odd: " + std::to_string(lines.size()));

    if (lines.size() < 2)
        exit(std::to_string(lines.size()) + " line !!!");

    return lines;
}

static inline Sets lines_to_sets(const std::vector<std::string> &lines)
{
    Sets sets(1);

    auto current_id = std::stoi(std::string(lines[0].begin() + 4, lines[0].end()));
    for (size_t i = 0; i < lines.size(); i += 2)
    {
        auto id = std::stoi(std::string(lines[i].begin() + 4, lines[i].end()));
        if (id != current_id)
        {
            current_id = id;
            sets.resize(sets.size() + 1);
        }
        sets.back().push_back(lines[i + 1]);
    }

    return sets;
}

static inline Set lines_to_sequences(const std::vector<std::string> &lines)
{
    Set set(lines.size() / 2);

    for (size_t i = 0; i < lines.size(); i += 2)
        set[i / 2] = lines[i + 1];

    return set;
}

std::vector<std::string> read_raw_sequence_file(const std::filesystem::path &filename)
{
    MappedFile mappedFile(filename);
    return split_lines(mappedFile.getData(), mappedFile.getSize());
}

Sets read_set_fasta(const std::filesystem::path &filename)
{
    std::vector<std::string> lines = read_raw_sequence_file(filename);
    check_lines_number(lines);
    return lines_to_sets(lines);
}

Set read_seq_fasta(const std::filesystem::path &filename)
{
    std::vector<std::string> lines = read_raw_sequence_file(filename);
    check_lines_number(lines);
    return lines_to_sequences(lines);
}