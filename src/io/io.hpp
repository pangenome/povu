#ifndef IO_HPP
#define IO_HPP

#include <cstddef>
#include <string>
#include <vector>

#include "../cli/app.hpp"
#include "../graph/bidirected.hpp"
#include "../graph/tree.hpp"
#include "../genomics/genomics.hpp"

#include "../graph/graph.hpp"


namespace io::from_gfa {
namespace bd = povu::bidirected;

povu::graph::Graph to_pv_graph(const char *filename, const core::config& app_config);
bd::VariationGraph to_bd(const char* filename, const core::config& app_config);
}; // namespace io::from_gfa

namespace povu::io::generic {
namespace pgt = povu::graph_types;

#define FILE_ERROR(name)                                                       \
  {                                                                            \
    std::string e = "Error, Failed to open the file " + name; \
    throw std::invalid_argument(e);                                              \
  }

std::vector<std::filesystem::path> get_files(const std::string& dir_path, const std::string& ext);

void write_txt(const std::vector<pgt::flubble>& flubbles, const std::string& base_name, const core::config& app_config);
}; // namespace io::generic


namespace povu::io::bub {
using povu::graph_types::id_n_orientation_t;
using povu::graph_types::id_n_cls;
namespace pvtr = povu::tree;
namespace pgt = povu::graph_types;

void write_bub(const pvtr::Tree<pgt::flubble>& bt, const std::string& base_name, const core::config& app_config);
/**
  * @brief Read a flb file but only return the canonical flubbles
 */
std::vector<pgt::flubble> read_canonical_fl(const std::string& fp);
} // namespace povu::io::bub

namespace povu::io::vcf {
using povu::genomics::vcf::vcf_record;



void write_vcfs(const std::map<std::size_t,
                std::vector<vcf_record>>& vcf_records,
                const bidirected::VariationGraph& bd_vg,
                const core::config& app_config);
} // namespace io::vcf


#endif
