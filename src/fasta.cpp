/*
 * Copyright 2022 - UPMEM
 */

#include <algorithm>
#include <fstream>

#include "fasta.hpp"

class File
{
  FILE *f = nullptr;

public:
  File() = delete;

  /// @brief Open a file
  /// @param filename Name of file to open
  File(const std::string &filename) : f(fopen(filename.c_str(), "r")) {}
  File(const File &) = delete;
  File(File &&) = delete;
  File &operator=(const File &) = delete;
  File &operator=(File &&) = delete;

  /// @brief Return true if file is valid
  inline operator bool() { return f != nullptr; }

  /// @brief Jump at the end of the file
  inline void SeekEnd() { fseek(f, 0, SEEK_END); }

  /// @brief Jump to the beginning of the file
  inline void Rewind() { ::rewind(f); }

  /// @brief Returns current position in the file
  /// @return
  inline size_t Tell() { return ftell(f); }

  /// @brief Reads the file
  /// @param ptr Pointer to a buffer the data will be writen to
  /// @param size sizeof data
  /// @param n Number of data to read
  /// @return Number of data read
  inline size_t Read(auto &ptr, size_t n)
  {
    return fread(ptr.data(), sizeof(ptr[0]), n, f);
  }

  ~File()
  {
    if (*this)
      fclose(f);
  }
};

std::string file_to_memory(const std::filesystem::path &filename)
{
  File file(filename.c_str());
  if (!file)
    exit("Invalid filename: " + filename.native());

  file.SeekEnd();
  size_t size = file.Tell();

  std::string buffer(size, '\0');

  file.Rewind();

  if (file.Read(buffer, size) != size)
    exit("Error reading: " + filename.native());

  return buffer;
}

std::vector<std::string> split_lines(const std::string &mapped_file)
{
  auto line_count = std::count(mapped_file.begin(), mapped_file.end(), '\n');

  if (mapped_file.back() != '\n')
    line_count++;

  std::vector<std::string> lines(line_count);
  size_t line = 0;
  size_t line_begin = 0;
  size_t line_size = 0;
  for (const auto &c : mapped_file)
  {
    if (c != '\n')
      ++line_size;
    else
    {
      if (line_size == 0)
        break;
      lines[line].resize(line_size);
      std::copy(
          &mapped_file.at(line_begin),
          &mapped_file.at(line_begin + line_size),
          lines[line].begin());
      line++;
      line_begin += line_size + 1;
      line_size = 0;
    }
  }

  if (line_size != 0)
  {
    lines[line].resize(line_size);
    std::copy(
        &mapped_file.at(line_begin),
        &mapped_file[line_begin + line_size],
        lines[line].begin());
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
  return file_to_memory(filename) |
         split_lines;
}

Sets read_set_fasta(const std::filesystem::path &filename)
{
  return read_raw_sequence_file(filename) |
         check_lines_number |
         lines_to_sets;
}

Set read_seq_fasta(const std::filesystem::path &filename)
{
  return read_raw_sequence_file(filename) |
         check_lines_number |
         lines_to_sequences;
}
