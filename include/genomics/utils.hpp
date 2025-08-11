#ifndef POVU_GRAPH_UTILS_HPP
#define POVU_GRAPH_UTILS_HPP

#include <algorithm>
#include <deque>
#include <format>
#include <set>
#include <stack>
#include <string>
#include <sys/types.h>
#include <utility>
#include <vector>

#include "../graph/bidirected.hpp"
#include "../common/types/genomics.hpp"
#include "../common/types/pvst.hpp"
#include "../common/types/core.hpp"
#include "../common/types/graph.hpp"


namespace povu::genomics::utils {
inline constexpr std::string_view MODULE = "povu::graph_utils";

using namespace povu::types::graph;
namespace pt = povu::types;
namespace bd = povu::bidirected;
namespace pvt = povu::types::genomics;

namespace pvst = povu::types::pvst;
namespace pc = povu::constants;
namespace pgt = povu::types::graph;




namespace variants {
  void comp_itineraries(const bd::VG &g, const pvt::walk &w, pt::idx_t w_idx, pvt::Exp &rw);
} // namespace variants

namespace graph {
void find_walks(const bd::VG &g, pvt::RoV &rov);
} // namespace graph

} // namespace povu::genomics::utils

#endif
