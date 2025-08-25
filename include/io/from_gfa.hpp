#ifndef PV_IO_FRM_GFA_HPP
#define PV_IO_FRM_GFA_HPP

#include <fstream> // for std::ifstream
#include <gfa.h> // from liteseq
#include <ostream>
#include <string>
#include <vector>

#include "../../include/common/compat.hpp"
#include "../../app/cli/app.hpp"
#include "../../include/common/utils.hpp"
#include "../../include/graph/bidirected.hpp"


namespace povu::io::from_gfa {
inline constexpr std::string_view MODULE = "povu::io::from_gfa";

namespace lq = liteseq;
namespace bd = povu::bidirected;
namespace pgt = povu::types::graph;
namespace pc = povu::constants;

bd::VG *to_bd(const core::config& app_config);
}; // namespace io::from_gfa
#endif
