#ifndef PVST_IO_HPP
#define PVST_IO_HPP

#include <cstddef>
#include <format>
#include <fstream> // for std::ifstream
#include <gfa.h>   // from liteseq
#include <ostream>
#include <string>
#include <vector>


//#include "../../include/common/utils.hpp"
//#include "../../include/graph/bidirected.hpp"
#include "../../include/graph/tree.hpp"
#include "../cli/app.hpp"
#include "./common.hpp"

namespace povu::io::pvst {
using povu::graph_types::id_n_cls;
using povu::graph_types::id_or_t;
namespace pvtr = povu::tree;
namespace pgt = povu::graph_types;
namespace pvst = povu::types::pvst;
namespace pc = povu::constants;
namespace pu = povu::utils;

void write_bub(const pvtr::Tree<pvst::Vertex> &bt, const std::string &base_name,
               const core::config &app_config);
/**
 * @brief Read a flb file but only return the canonical flubbles
 */
std::vector<pgt::flubble> read_canonical_fl(const std::string &fp);
} // namespace povu::io::pvst

#endif
