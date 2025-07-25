#ifndef POVU_GRAPH_UTILS_HPP
#define POVU_GRAPH_UTILS_HPP

#include <algorithm>
#include <format>
#include <vector>

#include "./types/pvst.hpp"
#include "./types/genomics.hpp"
#include "../graph/bidirected.hpp"
#include "types/graph.hpp"


namespace povu::graph_utils {

inline constexpr std::string_view MODULE = "povu::graph_utils";
using namespace povu::types::graph;
namespace pc = povu::constants;
namespace pt = povu::types;
namespace pgt = povu::types::graph;
namespace bd = povu::bidirected;
namespace pvst = povu::types::pvst;
namespace pvt = povu::types::genomics;

// Maximum number of steps to take from flubble start to end
const pt::idx_t MAX_FLUBBLE_STEPS{20};

void get_walks(const bd::VG &g, const pvst::VertexBase *pvst_vtx_ptr,
               std::vector<pgt::Walk>);

} // namespace povu::graph_utils

#endif
