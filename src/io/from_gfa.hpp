#ifndef PV_IO_FRM_GFA_HPP
#define PV_IO_FRM_GFA_HPP

#include <format>
#include <fstream> // for std::ifstream
#include <gfa.h>   // from liteseq
#include <ostream>
#include <string>
#include <vector>

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


#endif
