#ifndef PVST_IO_HPP
#define PVST_IO_HPP

#include <cstddef>
#include <fstream> // for std::ifstream
#include <gfa.h>
#include <ostream>
#include <string>
#include <vector>

#include "../../app/cli/app.hpp"
#include "../../include/common/types/compat.hpp"
#include "../../include/graph/tree.hpp"
#include "./common.hpp"

namespace povu::io::pvst {
using povu::types::graph::id_n_cls;
using povu::types::graph::id_or_t;
namespace pvtr = povu::tree;
namespace pgt = povu::types::graph;
namespace pvst = povu::types::pvst;
namespace pc = povu::constants;
namespace pu = povu::utils;
namespace pt = povu::types;

pvtr::Tree read_pvst(const std::string &fp);

void write_pvst(const pvtr::Tree &bt, const std::string &base_name, const core::config &app_config);
/**
 * @brief Read a flb file but only return the canonical flubbles
 */
std::vector<pgt::flubble> read_canonical_fl(const std::string &fp);
} // namespace povu::io::pvst

#endif
