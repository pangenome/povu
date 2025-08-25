#ifndef IO_HPP
#define IO_HPP

#include <cstddef>
#include <fstream> // for std::ifstream
#include <gfa.h>   // from liteseq
#include <ostream>
#include <string>
#include <vector>

#include "../../app/cli/app.hpp"
#include "../../include/common/compat.hpp"
#include "../../include/common/utils.hpp"
#include "../../include/graph/bidirected.hpp"

namespace povu::io::common {
namespace pgt = povu::types::graph;
namespace fs = std::filesystem;

#define FILE_ERROR(name)                                                       \
  {                                                                            \
    std::string e = "Error, Failed to open the file " + name; \
    throw std::invalid_argument(e);                                              \
  }
/**
 * @brief Get the list of files in a dir with a given name
 *
 */
std::vector<fs::path> get_files(const std::string& dir_path, const std::string& ext);

void write_txt(const std::vector<pgt::flubble>& flubbles, const std::string& base_name, const core::config& app_config);

void read_lines_to_vec_str(const std::string &fp, std::vector<std::string> *v);

}; // namespace io::generic



#endif
