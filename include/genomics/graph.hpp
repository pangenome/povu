#ifndef POVU_GRAPH_UTILS_HPP
#define POVU_GRAPH_UTILS_HPP

#include <algorithm>
#include <deque>
#include <set>
#include <stack>
#include <string>
#include <sys/types.h>
#include <utility>
#include <vector>

#include "../common/compat.hpp"
#include "../common/types/core.hpp"
#include "../common/types/genomics.hpp"
#include "../common/types/graph.hpp"
#include "../common/types/pvst.hpp"
#include "../common/log.hpp"
#include "../graph/bidirected.hpp"

namespace povu::genomics::graph {
inline constexpr std::string_view MODULE = "povu::genomics::graph";

namespace pvt = povu::types::genomics;
namespace pvst = povu::types::pvst;
namespace pgt = povu::types::graph;

void find_walks(const bd::VG &g, pvt::RoV &rov);
} // namespace povu::genomics::graph

#endif
