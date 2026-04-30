#include "ita/variation/color.hpp"

#include <optional> // for optional, operator==

#include <quilt/graph_types.hpp> // for v_end_e, side_n_id_t, side_n_idx_t
#include <quilt/shim.hpp>	 // for contains
#include <quilt/types.hpp>	 // for qt

#include "povu/graph/bidirected.hpp" // for VG, bd
#include "povu/graph/pvst.hpp"	     // for Tree, VertexBase

namespace ita::color
{
namespace pvst = oza::pvst;

bool has_any_refs(const bd::VG &g, const std::set<qt::id_t> &to_call_ref_ids,
		  qt::u32 s_v_idx, qt::u32 t_v_idx)
{
	for (qt::u32 r_idx : to_call_ref_ids) {
		const std::vector<qt::idx_t> &s_idxs =
			g.get_vertex_ref_idxs(s_v_idx, r_idx);

		const std::vector<qt::idx_t> &t_idxs =
			g.get_vertex_ref_idxs(t_v_idx, r_idx);

		if (!s_idxs.empty() || !t_idxs.empty())
			return true;
	}

	return false;
}

/**
 * compute a DFS from the root to find vertices to color
 */
std::set<qt::u32> color_pvst(const bd::VG &g, const pvst::Tree &pvst,
			     const std::set<qt::id_t> &to_call_ref_ids)
{
	std::stack<qt::u32> s;
	s.push(pvst.root_idx());

	std::set<qt::u32> to_call;

	while (!s.empty()) {
		qt::u32 curr_pvst_v_idx = s.top();
		s.pop();

		const pvst::VertexBase *curr_v_ptr =
			pvst.get_vertex_const_ptr(curr_pvst_v_idx);

		if (pvst::to_clan(curr_v_ptr->get_fam()) ==
		    pvst::vc_e::subflubble)
			continue;

		// if the current vertex has route params
		if (curr_v_ptr->get_route_params() != std::nullopt) {
			auto [l, r, _] = curr_v_ptr->get_route_params().value();
			auto [start_id, __] = l;
			auto [stop_id, ___] = r;

			qt::u32 s_v_idx = g.v_id_to_idx(start_id);
			qt::u32 t_v_idx = g.v_id_to_idx(stop_id);

			// color this vertex if either the start or stop vertex
			// has any of the refs we want to call
			if (has_any_refs(g, to_call_ref_ids, s_v_idx, t_v_idx))
				to_call.insert(curr_pvst_v_idx);

			qt::u32 p_pvst_v_idx =
				pvst.get_parent_idx(curr_pvst_v_idx);

			// remove parent from to_call if both parent and current
			// vertex are in to_call.
			if (qs::contains(to_call, p_pvst_v_idx) &&
			    qs::contains(to_call, curr_pvst_v_idx))
				to_call.erase(p_pvst_v_idx);
		}

		// push children to stack
		const std::vector<qt::idx_t> &children_idxs =
			pvst.get_children(curr_pvst_v_idx);

		for (qt::idx_t c_idx : children_idxs)
			s.push(c_idx);
	}

	return to_call;

	// std::set<qt::u32> colored_vtxs;
	// return colored_vtxs;
}

} // namespace ita::color
