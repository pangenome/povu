#ifndef FROM_PVST_IO_HPP
#define FROM_PVST_IO_HPP

#include <cstddef>
#include <fstream> // for std::ifstream
#include <gfa.h>
#include <ostream>
#include <string>
#include <vector>


#include <fmt/color.h>

#include "../../app/cli/app.hpp"
#include "../../include/common/compat.hpp"
#include "../../include/graph/tree.hpp"
#include "./common.hpp"
#include "../../include/common/log.hpp"

namespace povu::io::from_pvst {
constexpr std::string_view MODULE = "povu::io::from_pvst";

using povu::types::graph::id_n_cls;
using povu::types::graph::id_or_t;
namespace pvtr = povu::tree;
namespace pgt = povu::types::graph;
namespace pvst = povu::types::pvst;
namespace pc = povu::constants;
namespace pu = povu::utils;

pvtr::Tree read_pvst(const std::string &fp);


} // namespace povu::io::pvst

// add namespace alias for povu::compat
// below tells clang-tidy to skip that specific check for the next line.
// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pv_frm_pvst = povu::io::from_pvst;


#endif
