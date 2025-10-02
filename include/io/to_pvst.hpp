#ifndef TO_PVST_IO_HPP
#define TO_PVST_IO_HPP

#include <cstddef>
#include <fstream> // for std::ifstream
#include <liteseq/gfa.h>
#include <ostream>
#include <string>
#include <vector>

#include "../../app/cli/app.hpp"
#include "../common/compat.hpp"
#include "../graph/pvst.hpp"
#include "./common.hpp"

namespace povu::io::to_pvst {
using povu::types::graph::id_n_cls;
using povu::types::graph::id_or_t;
namespace pgt = povu::types::graph;
namespace pvst = povu::pvst;
namespace pc = povu::constants;
namespace pu = povu::utils;

void write_pvst(const pvst::Tree &bt, const std::string &base_name, const core::config &app_config);

} // namespace povu::io::pvst

// add namespace alias for povu::compat
// below tells clang-tidy to skip that specific check for the next line.
// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pv_to_pvst = povu::io::to_pvst;

#endif
