// TODO: remove this module

#ifndef POVU_GENOMICS_GRAPH_HPP
#define POVU_GENOMICS_GRAPH_HPP

#include <string_view> // for string_view

#include <oza/graph/bidirected.hpp> // for VG, bd
#include <oza/graph/pvst.hpp>	    // for VertexBase
#include <quilt/graph_types.hpp>    // for v_end_e, side_n_id_t, side_n_idx_t
#include <quilt/types.hpp>	    // for qt

#include "ita/variation/rov.hpp" // for RoV

namespace povu::genomics::graph
{
namespace pvst = oza::pvst;
namespace pgt = quilt::types::graph;

// Maximum number of steps to take from flubble start to end
const qt::u32 MAX_FLUBBLE_STEPS{1000};
const qt::u32 MAX_UNBLOCK_CTR{10000};

// direction for traversing a vertex in a bidirected graph
enum class dir_e : qt::u8 {
	in,
	out
};

std::string_view to_str(dir_e d);
const dir_e IN = dir_e::in;
const dir_e OUT = dir_e::out;
auto format_as(dir_e d);

typedef pgt::id_or_t idx_or_t; // specifically for idx instead of id

qt::status_t find_walks(const bd::VG &g, ir::RoV &rov);
} // namespace povu::genomics::graph

#endif // POVU_GENOMICS_GRAPH_HPP
