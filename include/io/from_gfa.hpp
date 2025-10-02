#ifndef PV_IO_FRM_GFA_HPP
#define PV_IO_FRM_GFA_HPP

#include "../../app/cli/app.hpp"
#include "../../include/common/compat.hpp"
#include "../../include/common/utils.hpp"
#include "../../include/graph/bidirected.hpp"
#include <fstream> // for std::ifstream
#include <liteseq/gfa.h>   // from liteseq
#include "../common/progress.hpp"

#include <ostream>
#include <string>
#include <vector>

namespace povu::io::from_gfa {
inline constexpr std::string_view MODULE = "povu::io::from_gfa";

using namespace povu::progress;
namespace pv_prog = povu::progress;

namespace lq = liteseq;
namespace bd = povu::bidirected;
namespace pgt = povu::types::graph;
namespace pc = povu::constants;

/**
 * Remember to free the returned graph after use
 */
bd::VG *to_bd(const core::config& app_config);
}; // namespace io::from_gfa
#endif
