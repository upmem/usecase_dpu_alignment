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
  File(const std::string &filename) : f(fopen(filename.c_str(), "r")) {}
  File(const File &) = delete;
  File(File &&) = delete;
  File &operator=(const File &) = delete;
  File &operator=(File &&) = delete;

  inline bool is_valid() { return f != nullptr; }
  inline bool is_invalid() { return f == nullptr; }
  inline void seek_end() { fseek(f, 0, SEEK_END); }
  inline void rewind() { ::rewind(f); }
  inline size_t tell() { return ftell(f); }
  inline size_t read(void *__restrict__ ptr, size_t size, size_t n)
  {
    return fread(ptr, size, n, f);
  }

  ~File()
  {
    if (is_valid())
      fclose(f);
  }
};

std::string file_to_memory(const std::filesystem::path &filename)
{
  File f(filename.c_str());
  if (f.is_invalid())
    exit("Invalid filename: " + filename.native());

  f.seek_end();
  size_t size = f.tell();

  std::string buffer(size, '\0');

  f.rewind();

  if (f.read(buffer.data(), sizeof(char), size) != size)
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

static inline Set lines_to_sequences(const Sequences &lines)
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