#include "povu/variation/rov.hpp"

#include <cstddef>  // for size_t
#include <cstdlib>  // for exit, EXIT_FAILURE
#include <iterator> // for back_insert_iterator, bac...
#include <map>	    // for map
#include <optional> // for optional, operator==
#include <string>   // for basic_string, string
#include <utility>  // for move
#include <vector>   // for vector

// #include "fmt/core.h"		  // for format_to
// #include "indicators/setting.hpp" // for PostfixText
// #include "povu/common/app.hpp"	  // for config
#include "povu/common/constants.hpp"
#include "povu/common/core.hpp" // for pt
// #include "povu/common/log.hpp"
// #include "povu/common/progress.hpp" // for set_progress_bar_com...
// #include "povu/common/utils.hpp"
#include "povu/genomics/graph.hpp" // for RoV, find_walks, pgt
// #include "povu/genomics/vcf.hpp"
#include "povu/graph/pvst.hpp" // for Tree, VertexBase
#include "povu/graph/types.hpp"
#include "povu/variation/color.hpp"

// #include "povu/common/utils.hpp" // for print_with_comma

namespace povu::var::rov
{
namespace pgg = povu::genomics::graph;
namespace pvst = povu::pvst;
// using namespace povu::progress;

constexpr var_type_e ins = var_type_e::ins;
constexpr var_type_e del = var_type_e::del;
constexpr var_type_e sub = var_type_e::sub;

// ------------

std::ostream &operator<<(std::ostream &os, var_type_e vt)
{
	return os << to_string_view(vt);
}

var_type_e covariant(var_type_e a) noexcept
{
	switch (a) {
	case var_type_e::ins:
		return var_type_e::del;
	case var_type_e::del:
		return var_type_e::ins;
	case var_type_e::sub:
		return var_type_e::sub;
	default:
		return var_type_e::und;
	}
}

bool is_valid(pt::u32 i)
{
	return i != pc::INVALID_IDX;
}

/**
 * map each step to its position in the walk
 */
std::map<ptg::step_t, pt::u32> positions(const ptg::walk_t &w)
{
	std::map<ptg::step_t, pt::u32> m;
	for (pt::u32 i{}; i < w.size(); i++)
		m[w[i]] = i;

	return m;
}

/**
 * look up each step in the walk in the map
 */
std::vector<pt::u32> comp_lookup(const ptg::walk_t &w,
				 std::map<ptg::step_t, pt::u32> &m)
{
	std::vector<pt::u32> lu(w.size(), pc::INVALID_IDX);
	for (pt::u32 i{}; i < w.size(); i++)
		if (pv_cmp::contains(m, w[i]))
			lu[i] = m[w[i]];

	return lu;
}

pt::slice_t find_context(const std::vector<pt::u32> &w, pt::u32 i)
{
	pt::u32 left = i;
	pt::u32 right = i + 1;

	while (right < w.size() && !is_valid(w[right]))
		right++;

	return {left, right - left};
}

void find_rovs(const std::vector<pt::u32> &lu, pairwise_variants &pv)
{
	auto is_ins = [&](pt::u32 i, const pt::slice_t &sl) -> bool
	{
		return (i > 0) && (lu[i + sl.len()] - lu[i - 1] == 1);
	};

	auto is_del = [&](pt::u32 i) -> bool
	{
		return i > 0 && is_valid(lu[i - 1]) && is_valid(lu[i]) &&
		       lu[i] - lu[i - 1] != 1;
	};

	auto find_alt_start = [&](const pt::slice_t &sl, var_type_e t,
				  pt::u32 i) -> pt::slice_t
	{
		pt::u32 alt_len = t == ins ? 0 : lu[i + sl.len()] - lu[i - 1];
		pt::u32 alt_start = lu[i - 1];

		if (t == sub || alt_start == 0) {
			if (t != ins)
				alt_len--;
			if (t != ins)
				alt_start++;
		}

		return {alt_start, alt_len};
	};

	// because we move left to right we will always
	// start with the leftmost side of the slice
	for (pt::u32 i{}; i < lu.size();) {
		if (!is_valid(lu[i])) { // sub or ins
			pt::slice_t sl = find_context(lu, i);
			var_type_e t = is_ins(i, sl) ? ins : sub;
			pt::slice_t alt_sl = find_alt_start(sl, t, i);
			pv.add_variant({sl, alt_sl, covariant(t)});
			i += sl.len();
			continue;
		}

		// del
		if (is_del(i)) {
			pt::slice_t sl{i - 1, 0};
			pt::slice_t alt_sl = {lu[i - 1] + 1,
					      lu[i] - lu[i - 1] - 1};
			pv.add_variant({sl, alt_sl, covariant(del)});
		}

		i++;
	}
}

pairwise_variants lineup_pairs(const ptg::walk_t &w1, const ptg::walk_t &w2,
			       pairwise_variants &pv)
{
	// std::map<ptg::step_t, pt::u32> pos_map1 = positions(w1);
	std::map<ptg::step_t, pt::u32> pos_map2 = positions(w2);

	std::vector<pt::u32> lu1 = comp_lookup(w1, pos_map2);
	// std::vector<pt::u32> lu2 = comp_lookup(w2, pos_map1);

	find_rovs(lu1, pv);
	// find_rovs(lu2, true, pv);

	return pv;
}

pairwise_variants compare_pair(const RoV &r, pt::u32 i, pt::u32 j)
{
	pairwise_variants pv(i, j);
	lineup_pairs(r.get_walk(i), r.get_walk(j), pv);

	return pv;
}

void find_hidden(RoV &r)
{
	// std::cerr << "RoV: " << r.as_str() << "\n";
	// for (int i = 0; i < r.get_walks().size(); i++) {
	//	std::cerr << i << "\t" << pgt::to_string(r.get_walk(i)) << "\n";
	// }

	//
	//  bool dbg = r.as_str() == ">1546>1551" ? true : false;
	//  dbg = false;
	//  for (int i = 0; i < r.get_walks().size(); i++) {
	//	if (dbg) {
	//		const ptg::walk_t &w = r.get_walk(i);
	//		std::cerr << i << ": " << pgt::to_string(w) << "\n";
	//	}
	//  }

	// const pt::u32 WC = r.walk_count();
	// std::set<pt::up_t<pt::u32>> seen;
	// for (pt::u32 i{}; i < WC; i++) {
	//	for (pt::u32 j{}; j < WC; j++) {
	//		auto p = pt::up_t<pt::u32>{i, j};

	//		if (i == j || pv_cmp::contains(seen, p))
	//			continue; // skip self or already seen

	//		r.add_irreducible(compare_pair(r, i, j));
	//		seen.insert(p);
	//	}
	// }

	std::set<pt::up_t<pt::u32>> flanks = find_non_planar(r.get_walks());
	r.extend_flanks(flanks);

	// if (dbg)
	//	for (auto &pv : r.get_irreducibles())
	//		std::cerr << pv << "\n";
}

// ------------

/**
 * Check if a vertex in the pvst is a flubble leaf
 * A flubble leaf is a vertex that has no children that are also
 * flubbles
 */
bool is_fl_leaf(const pvst::Tree &pvst, pt::idx_t pvst_v_idx) noexcept
{
	const pvst::VertexBase *pvst_v_ptr =
		pvst.get_vertex_const_ptr(pvst_v_idx);

	// we assume that the vertex has a clan
	pvst::vf_e prt_fam = pvst_v_ptr->get_fam();
	if (pvst::to_clan(prt_fam).value() != pvst::vc_e::fl_like) {
		return false; // not a flubble
	}

	for (pt::idx_t v_idx : pvst.get_children(pvst_v_idx)) {
		pvst::vf_e c_fam = pvst.get_vertex_const_ptr(v_idx)->get_fam();
		if (pvst::to_clan(c_fam) == pvst::vc_e::fl_like)
			return false;
	}

	return true;
}

/**
 * true when the vertex is a flubble leaf or a leaf in the pvst
 */
bool should_call(const pvst::Tree &pvst, const pvst::VertexBase *pvst_v_ptr,
		 pt::idx_t pvst_v_idx)
{
	if (auto opt_route_params = pvst_v_ptr->get_route_params())
		return is_fl_leaf(pvst, pvst_v_idx) || pvst.is_leaf(pvst_v_idx);

	return false;
}

void find_pvst_rovs(const bd::VG &g, const pvst::Tree &pvst,
		    const std::set<pt::u32> &colored_vtxs, std::vector<RoV> &rs)
{
	// auto can_be_non_planar = [](const pvst::vf_e fam) -> bool
	// {
	//	return fam == pvst::vf_e::flubble;
	// };

	for (pt::u32 i : colored_vtxs) { // i is pvst_v_idx
		// std::cout << "Colored pvst vertex idx: " << v_idx << "\n";
		const pvst::VertexBase *v = pvst.get_vertex_const_ptr(i);
		RoV r{v};

		// if (r.as_str() == ">2597>2621") {
		//	std::cerr << "Found target RoV\n";
		//	std::exit(EXIT_FAILURE);
		// }

		// get the set of walks for the RoV
		pt::status_t s = povu::genomics::graph::find_walks(g, r);

		// no walks found, skip this RoV
		if (r.get_walks().size() == 0 || s != 0)
			continue;

		if (r.can_be_non_planar())
			find_hidden(r);

		// if (r.get_pvst_vtx()->get_fam() != pvst::vf_e::tiny ||
		//     r.get_pvst_vtx()->get_fam() != pvst::vf_e::parallel)
		//	find_hidden(r);

		rs.emplace_back(std::move(r));
	}
}

// TODO: remove this function (unused)
void eval_vertex(const bd::VG &g, const pvst::Tree &pvst, pt::u32 pvst_v_idx,
		 std::vector<RoV> &rs)
{
	// std::string s = ">3597>3600";
	const pvst::VertexBase *pvst_v_ptr =
		pvst.get_vertex_const_ptr(pvst_v_idx);

	std::string ss = pvst_v_ptr->as_str();

	// if (s == ss) {
	//	std::cerr << "Evaluating pvst vertex: " << ss << "\n";
	// }

	if (should_call(pvst, pvst_v_ptr, pvst_v_idx)) {
		RoV r{pvst_v_ptr};

		// get the set of walks for the RoV
		pt::status_t s = povu::genomics::graph::find_walks(g, r);

		// no walks found, skip this RoV
		if (r.get_walks().size() == 0 || s != 0)
			return;

		find_hidden(r);
		rs.emplace_back(std::move(r));
	}
}

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

/**
 * find walks in the graph based on the leaves of the pvst
 * initialize RoVs from flubbles
 */
std::vector<RoV> gen_rov(const std::vector<pvst::Tree> &pvsts, const bd::VG &g,
			 const std::set<pt::id_t> &to_call_ref_ids)
{
	// the set of RoVs to return
	std::vector<RoV> rs;
	rs.reserve(pvsts.size());

	// // reuse msg buffer for progress bar
	// DynamicProgress<ProgressBar> bars;
	// std::string prog_msg;
	// prog_msg.reserve(128);

	// auto update_prog = [&](pt::idx_t i, pt::idx_t v_idx, pt::idx_t total,
	//		       std::size_t bar_idx)
	// {
	//	prog_msg.clear();
	//	fmt::format_to(std::back_inserter(prog_msg),
	//		       "Generating RoVs for PVST {} ({}/{})", i + 1,
	//		       v_idx + 1, total);
	//	bars[bar_idx].set_option(option::PostfixText{prog_msg});
	//	bars[bar_idx].set_progress(v_idx + 1);
	// };

	for (pt::idx_t i{}; i < pvsts.size(); i++) { // for each pvst
		const pvst::Tree &pvst = pvsts.at(i);
		// loop through each tree
		// const pt::idx_t total = pvst.vtx_count();

		std::set<pt::u32> colored_vtxs =
			pvc::color_pvst(g, pvst, to_call_ref_ids);

		// std::cerr << "clrd " << colored_vtxs.size() << "\n";

		find_pvst_rovs(g, pvst, colored_vtxs, rs);
		// reset progress bar
		// ProgressBar bar;
		// set_progress_bar_common_opts(&bar);
		// std::size_t bar_idx = bars.push_back(bar);
		// set_progress_bar_common_opts(&bar, pvst.vtx_count());

		// for (pt::u32 v_idx{}; v_idx < pvst.vtx_count(); v_idx++) {
		//	// update progress bar
		//	// if (app_config.show_progress())
		//	//	update_prog(i, v_idx, total, bar_idx);

		//	eval_vertex(g, pvst, v_idx, rs);

		//	if (app_config.show_progress()) //
		//	//	bars[bar_idx].mark_as_completed();
		// }
	}

	return rs;
}

} // namespace povu::var::rov
