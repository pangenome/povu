#ifndef IO_HPP
#define IO_HPP

#include <cstddef>
#include <string>
#include <vector>
#include <format>
#include <cstddef>
#include <ostream>
#include <fstream> // for std::ifstream
#include <string>
#include <vector>
#include <gfa.h> // from liteseq

#include "../cli/app.hpp"
#include "../../include/graph/bidirected.hpp"
#include "../../include/graph/tree.hpp"
#include "../../include/common/utils.hpp"
#include "../../include/graph/tree.hpp"


namespace povu::io::from_gfa {
namespace lq = liteseq;
namespace bd = povu::bidirected;
namespace pgt = povu::graph_types;
namespace pt = povu::types;

bd::VG *to_bd(const char* filename, const core::config& app_config);
}; // namespace io::from_gfa

namespace povu::io::common {
namespace pgt = povu::graph_types;
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

namespace povu::io::bub {
using povu::graph_types::id_or_t;
using povu::graph_types::id_n_cls;
namespace pvtr = povu::tree;
namespace pgt = povu::graph_types;

void write_bub(const pvtr::Tree<pgt::flubble>& bt, const std::string& base_name, const core::config& app_config);
/**
  * @brief Read a flb file but only return the canonical flubbles
 */
std::vector<pgt::flubble> read_canonical_fl(const std::string& fp);
} // namespace povu::io::bub

/*
namespace povu::io::vcf {
using povu::genomics::vcf::vcf_record;



void write_vcfs(const std::map<std::size_t,
                std::vector<vcf_record>>& vcf_records,
                const bidirected::VariationGraph& bd_vg,
                const core::config& app_config);
} // namespace io::vcf
*/

#endif
