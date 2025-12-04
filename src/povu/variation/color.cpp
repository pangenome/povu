#include "povu/variation/color.hpp"

#include <optional> // for optional, operator==
#include <vector>   // for vector

#include "povu/common/constants.hpp"
#include "povu/common/core.hpp"	     // for pt
#include "povu/graph/bidirected.hpp" // for VG, bd
#include "povu/graph/pvst.hpp"	     // for Tree, VertexBase
#include "povu/graph/types.hpp"

namespace povu::var::color
{
namespace pvst = povu::pvst;

// using namespace povu::progress;

bool has_all_refs(const bd::VG &g, const std::set<pt::id_t> &to_call_ref_ids,
		  pt::u32 s_v_idx, pt::u32 t_v_idx)
{
	for (pt::u32 r_idx : to_call_ref_ids) {
		const std::vector<pt::idx_t> &s_idxs =
			g.get_vertex_ref_idxs(s_v_idx, r_idx);

		const std::vector<pt::idx_t> &t_idxs =
			g.get_vertex_ref_idxs(t_v_idx, r_idx);

		if (s_idxs.empty() || t_idxs.empty())
			return false;
	}

	return true;
}

bool has_any_refs(const bd::VG &g, const std::set<pt::id_t> &to_call_ref_ids,
		  pt::u32 s_v_idx, pt::u32 t_v_idx)
{
	for (pt::u32 r_idx : to_call_ref_ids) {
		const std::vector<pt::idx_t> &s_idxs =
			g.get_vertex_ref_idxs(s_v_idx, r_idx);

		const std::vector<pt::idx_t> &t_idxs =
			g.get_vertex_ref_idxs(t_v_idx, r_idx);

		if (s_idxs.empty() || t_idxs.empty())
			return false;
	}

	return true;
}

// given a vertex go up the tree until you find the flubble leaf which
// has all to call ref ids
std::optional<pt::u32> find_vertex(const bd::VG &g, const pvst::Tree &pvst,
				   const std::set<pt::id_t> &to_call_ref_ids,
				   pt::u32 pvst_v_idx)
{
	// std::string s = ">3597>3600";
	// std::string ss = pvst.get_vertex(pvst_v_idx).as_str();
	// auto dbg = (ss == s) ? true : false;

	auto is_parent_valid = [&](pt::u32 parent_pvst_v_idx)
	{
		return parent_pvst_v_idx != pvst.root_idx() &&
		       parent_pvst_v_idx != pc::INVALID_IDX &&
		       // if the parent has too many children, skip
		       // avoids going too far up the tree
		       pvst.get_children(parent_pvst_v_idx).size() < 20;
	};

	do {
		const pvst::VertexBase *v =
			pvst.get_vertex_const_ptr(pvst_v_idx);

		if (v->get_route_params() == std::nullopt)
			return std::nullopt;

		auto [l, r, _] = v->get_route_params().value();
		auto [start_id, __] = l;
		auto [stop_id, ___] = r;

		pt::u32 s_v_idx = g.v_id_to_idx(start_id);
		pt::u32 t_v_idx = g.v_id_to_idx(stop_id);

		if (has_all_refs(g, to_call_ref_ids, s_v_idx, t_v_idx))
			return pvst_v_idx;

		pvst_v_idx = pvst.get_parent_idx(pvst_v_idx);
	} while (is_parent_valid(pvst_v_idx));

	return std::nullopt;
}

std::set<pt::u32> color_pvst_old(const bd::VG &g, const pvst::Tree &pvst,
				 const std::set<pt::id_t> &to_call_ref_ids)
{

	std::set<pt::u32> colored_vtxs;
	for (pt::u32 i{}; i < pvst.vtx_count(); i++) // i is pvst_v_idx
		if (pvst.is_leaf(i))
			if (auto j = find_vertex(g, pvst, to_call_ref_ids, i))
				colored_vtxs.insert(*j);

	return colored_vtxs;
}

std::set<pt::u32> color_pvst(const bd::VG &g, const pvst::Tree &pvst,
			     const std::set<pt::id_t> &to_call_ref_ids)
{
	// compute a DFS from the root to find vertices to color
	// pt::u32 r = pvst.root_idx();

	std::stack<pt::u32> s;
	s.push(pvst.root_idx());

	std::set<pt::u32> to_call;

	while (!s.empty()) {
		pt::u32 curr_pvst_v_idx = s.top();
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

			pt::u32 s_v_idx = g.v_id_to_idx(start_id);
			pt::u32 t_v_idx = g.v_id_to_idx(stop_id);

			if (has_all_refs(g, to_call_ref_ids, s_v_idx,
					 t_v_idx)) {
				// color this vertex
				to_call.insert(curr_pvst_v_idx);
			}

			pt::u32 p_pvst_v_idx =
				pvst.get_parent_idx(curr_pvst_v_idx);

			if (pv_cmp::contains(to_call, p_pvst_v_idx) &&
			    pv_cmp::contains(to_call, curr_pvst_v_idx)) {
				// remove parent from to_call
				to_call.erase(p_pvst_v_idx);
			}
		}

		// push children to stack
		const std::vector<pt::idx_t> &children_idxs =
			pvst.get_children(curr_pvst_v_idx);

		for (pt::idx_t c_idx : children_idxs)
			s.push(c_idx);
	}

	return to_call;

	// std::set<pt::u32> colored_vtxs;
	// return colored_vtxs;
}

} // namespace povu::var::color
