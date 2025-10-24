#ifndef POVU_GENOMICS_GRAPH_HPP
#define POVU_GENOMICS_GRAPH_HPP

#include <string>      // for basic_string, string
#include <string_view> // for string_view
#include <utility>     // for move
#include <vector>      // for vector

#include "povu/common/core.hpp"	     // for idx_t, pt
#include "povu/genomics/rov.hpp"     // for RoV
#include "povu/graph/bidirected.hpp" // for VG, bd
#include "povu/graph/pvst.hpp"	     // for VertexBase
#include "povu/graph/types.hpp"	     // for walk_t, id_or_t

namespace povu::genomics::graph
{
inline constexpr std::string_view MODULE = "povu::genomics::graph";
namespace pvst = povu::pvst;
namespace pgt = povu::types::graph;

// Maximum number of steps to take from flubble start to end
const pt::idx_t MAX_FLUBBLE_STEPS{20};

// direction for traversing a vertex in a bidirected graph
enum class dir_e : pt::u8 {
	in,
	out
};

std::string_view to_str(dir_e d);
const dir_e IN = dir_e::in;
const dir_e OUT = dir_e::out;
auto format_as(dir_e d);

typedef pgt::id_or_t idx_or_t; // specifically for idx instead of id

pt::status_t find_walks(const bd::VG &g, pgr::RoV &rov);
} // namespace povu::genomics::graph

#endif // POVU_GENOMICS_GRAPH_HPP
