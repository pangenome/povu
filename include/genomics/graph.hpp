#ifndef POVU_GENOMICS_GRAPH_HPP
#define POVU_GENOMICS_GRAPH_HPP

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

// Maximum number of steps to take from flubble start to end
const pt::idx_t MAX_FLUBBLE_STEPS{20};

// direction for traversing a vertex in a bidirected graph
enum class dir_e { in, out };

const dir_e IN = dir_e::in;
const dir_e OUT = dir_e::out;

typedef pgt::id_or_t idx_or_t; // specifically for idx instead of id

void find_walks(const bd::VG &g, pvt::RoV &rov);
} // namespace povu::genomics::graph

#endif // POVU_GENOMICS_GRAPH_HPP
