#include "povu/overlay/overlay.hpp"

#include <algorithm>
#include <csignal> // for raise, SIGINT
#include <cstdint>
#include <cstdio>
#include <cstdlib> // for exit, EXIT_FAILURE
#include <iostream>
#include <liteseq/refs.h> // for ref_walk, ref
#include <map>		  // for map
#include <string>
#include <utility>
#include <vector> // for vector

#include "povu/common/compat.hpp"
#include "povu/common/constants.hpp"
#include "povu/common/core.hpp"
#include "povu/common/log.hpp"
#include "povu/common/utils.hpp"

#include "povu/genomics/allele.hpp"

#include "povu/graph/types.hpp"
#include "povu/variation/rov.hpp"

namespace povu::overlay
{
namespace lq = liteseq;

constexpr pvr::var_type_e ins = pvr::var_type_e::ins;
constexpr pvr::var_type_e del = pvr::var_type_e::del;
constexpr pvr::var_type_e sub = pvr::var_type_e::sub;
constexpr pgt::or_e fo = pgt::or_e::forward;
constexpr pgt::or_e ro = pgt::or_e::reverse;

const pt::u8 SLICE_A_IDX{0};
const pt::u8 SLICE_B_IDX{1};
using prefix_sum = std::vector<pt::u32>;
using rs_to_ps = std::pair<pt::u32, prefix_sum>; // ref start to prefix sum

struct overlay_prefix_sum_t {
	pt::u32 ref_start_idx;
	prefix_sum p_sum_fwd;
	prefix_sum p_sum_rev;
};

class Overlays
{
	// i is the walk idx and j is the ref idx
	std::vector<std::vector<std::vector<overlay_prefix_sum_t>>>
		prefix_sums_;

	std::vector<sub_inv> sub_invs_;
	// key is walk idx and value is idx in overlays vector
	// std::map<pt::u32, pt::u32> walk_to_overlay_;

public:
	Overlays(pt::u32 walk_count, pt::u32 ref_count)
	    : prefix_sums_(walk_count,
			   std::vector<std::vector<overlay_prefix_sum_t>>(
				   ref_count)),
	      sub_invs_()
	{}

	[[nodiscard]]
	const std::vector<overlay_prefix_sum_t> &get_ps(pt::u32 w_idx,
							pt::u32 r_idx) const
	{
		return prefix_sums_.at(w_idx).at(r_idx);
	}

	void add_overlay(pt::u32 w_idx, pt::u32 r_idx,
			 std::vector<overlay_prefix_sum_t> &&p_sums)
	{
		prefix_sums_.at(w_idx).at(r_idx) = std::move(p_sums);
	}

	[[nodiscard]]
	std::vector<sub_inv> get_sub_invs_cpy() const
	{
		return sub_invs_;
	}

	[[nodiscard]]
	const std::vector<sub_inv> &get_sub_invs() const
	{
		return sub_invs_;
	}

	[[nodiscard]]
	bool has_sub_invs() const
	{
		return !sub_invs_.empty();
	}

	void add_sub_inv(sub_inv &&si)
	{
		sub_invs_.emplace_back(std::move(si));
	}
};

struct overlay_t {
	pt::idx_t graph_w_start_idx;
	pt::idx_t ref_start_idx;
	pt::idx_t len;
	ptg::or_e slice_or;
	pvr::var_type_e vt;
};

std::vector<pt::u32> comp_prefix_sum(const std::vector<pt::u8> &ov)
{
	std::vector<pt::u32> ps(ov.size(), 0);

	ps[0] = ov[0];
	for (pt::u32 i{1}; i < ov.size(); i++)
		ps[i] = ps[i - 1] + ov[i];

	// povu::utils::print_with_comma(std::cerr, ps, ',');
	// std::cerr << "\n";

	return ps;
}

std::vector<overlay_prefix_sum_t>
overlay_walk(const lq::ref_walk *ref_w, const pgt::walk_t &graph_w,
	     const std::vector<pt::idx_t> &ref_idx_starts)
{
	const pt::u32 GRAPH_W_LEN = graph_w.size();
	std::vector<pt::u8> ov(GRAPH_W_LEN, 0);
	std::vector<pt::u8> ov_rev(GRAPH_W_LEN, 0);

	// std::string w_str = pgt::to_string(graph_w);
	// std::cerr << w_str << "\t";
	// overlay_prefix_sum_t ops;
	// std::vector<rs_to_ps> prefix_sums;
	std::vector<overlay_prefix_sum_t> prefix_sums;
	prefix_sums.reserve(ref_idx_starts.size());

	for (pt::u32 ref_idx : ref_idx_starts) {
		pt::u32 ref_w_idx = ref_idx;
		pt::u32 step_idx{};
		for (; step_idx < graph_w.size(); step_idx++, ref_w_idx++) {
			auto [g_v_id, g_o] = graph_w[step_idx];

			pt::idx_t ref_v_id = ref_w->v_ids[ref_w_idx];
			pgt::or_e ref_o = ref_w->strands[ref_w_idx] ==
							  lq::strand::STRAND_FWD
						  ? pgt::or_e::forward
						  : pgt::or_e::reverse;

			if (ref_v_id != g_v_id || ref_o != g_o)
				ov[step_idx] = 1;
		}

		// std::cerr << "rw start " << ref_idx << "\n";
		// prefix_sums.emplace_back(ref_idx, comp_prefix_sum(ov));
		// reset ov for next ref start

		ref_w_idx = ref_idx;
		step_idx = 0;
		ref_w_idx++;
		for (; step_idx < graph_w.size() && ref_w_idx-- > 0;
		     step_idx++) {
			auto [g_v_id, g_o] = graph_w[step_idx];

			pt::idx_t ref_v_id = ref_w->v_ids[ref_w_idx];
			pgt::or_e ref_o = ref_w->strands[ref_w_idx] ==
							  lq::strand::STRAND_FWD
						  ? pgt::or_e::forward
						  : pgt::or_e::reverse;

			if (ref_v_id != g_v_id || ref_o != pgt::flip(g_o))
				ov_rev[step_idx] = 1;

			// if (ref_w_idx == 0)
			//	break;
		}

		// std::cerr << "rw start " << ref_idx << "\n";
		// prefix_sums.emplace_back(ref_idx, comp_prefix_sum(ov_rev));
		// reset ov for next ref start

		prefix_sums.push_back(overlay_prefix_sum_t{
			ref_idx, comp_prefix_sum(ov), comp_prefix_sum(ov_rev)});

		// std::cerr << "fw  ";
		// for (pt::u32 i{}; i < ov.size(); i++)
		//	std::cerr << static_cast<pt::u32>(ov[i]) << ", ";
		// std::cerr << "\n";

		// std::cerr << "ps  ";
		// auto x = comp_prefix_sum(ov);
		// for (pt::u32 i{}; i < ov.size(); i++)
		//	std::cerr << static_cast<pt::u32>(x[i]) << ", ";
		// std::cerr << "\n";

		// std::cerr << "rev ";
		// for (pt::u32 i{}; i < ov_rev.size(); i++)
		//	std::cerr << static_cast<pt::u32>(ov_rev[i]) << ", ";
		// std::cerr << "\n";

		// std::cerr << "ps rev  ";
		// x = comp_prefix_sum(ov_rev);
		// for (pt::u32 i{}; i < ov.size(); i++)
		//	std::cerr << static_cast<pt::u32>(x[i]) << ", ";
		// std::cerr << "\n";

		// reset ov for next ref start
		std::fill(ov.begin(), ov.end(), 0);
		std::fill(ov_rev.begin(), ov_rev.end(), 0);
	}

	// for (pt::u32 ref_idx : ref_idx_starts) {
	// }

	return prefix_sums;
}

/**
 * refs are in a walk if
 */
std::set<pt::id_t> find_refs_in_walk(const bd::VG &g, const pgt::walk_t &walk)
{
	pgt::step_t s = walk.front();
	pgt::step_t t = walk.back();

	auto [s_v_id, s_o] = s;
	auto [t_v_id, t_o] = t;

	const std::vector<std::vector<pt::idx_t>> &s_vtx_refs =
		g.get_vertex_refs(s_v_id);
	const std::vector<std::vector<pt::idx_t>> &t_vtx_refs =
		g.get_vertex_refs(t_v_id);

	std::set<pt::id_t> graph_walk_refs;

	for (pt::u32 ref_idx{}; ref_idx < g.ref_count(); ref_idx++) {
		if (s_vtx_refs[ref_idx].empty() || t_vtx_refs[ref_idx].empty())
			continue;

		graph_walk_refs.insert(ref_idx);
	}

	return graph_walk_refs;
}

std::map<pt::u32, std::set<pt::id_t>>
find_ref_walks(const bd::VG &g, const std::vector<pgt::walk_t> &walks)
{
	// map walks (walk idxs) to the refs that go through them
	std::map<pt::u32, std::set<pt::id_t>> refs_in_walks;
	for (pt::u32 w_idx{}; w_idx < walks.size(); ++w_idx)
		refs_in_walks[w_idx] = find_refs_in_walk(g, walks[w_idx]);

	return refs_in_walks;
}

std::optional<std::map<pgt::or_e, std::vector<pt::u32>>>
find_walk_sub_inv(const Overlays &overlays, const pt::u32 REF_COUNT,
		  const pt::u32 WALK_IDX)
{
	std::map<pgt::or_e, std::vector<pt::u32>> locs;

	for (pt::u32 r_idx{}; r_idx < REF_COUNT; ++r_idx) {
		const std::vector<overlay_prefix_sum_t> &pss =
			overlays.get_ps(WALK_IDX, r_idx);

		if (pss.size() != 1)
			continue;

		auto [_, ps_f, ps_r] = pss.front();

		if (ps_f.back() == 0) {
			locs[pgt::or_e::forward].push_back(r_idx);
		}
		else if (ps_r.back() == 0) {
			locs[pgt::or_e::reverse].push_back(r_idx);
		}
	}

	if (locs.size() != 2)
		return std::nullopt;

	return locs;
}

void find_sub_inv(Overlays &overlays, const pt::u32 REF_COUNT,
		  const pt::u32 WALK_COUNT)
{

	for (pt::u32 w_idx{}; w_idx < WALK_COUNT; ++w_idx) {

		auto locs_opt = find_walk_sub_inv(overlays, REF_COUNT, w_idx);

		if (!locs_opt.has_value())
			continue;

		// has both fwd and rev
		overlays.add_sub_inv(sub_inv{w_idx,
					     locs_opt->at(pgt::or_e::forward),
					     locs_opt->at(pgt::or_e::reverse)});
	}
}

/**
 * comp ps for each walk as maps to the ref
 */
Overlays overlay_walks(const bd::VG &g, const std::vector<pgt::walk_t> &walks)
{
	const pt::u32 REF_COUNT = g.ref_count();
	const pt::u32 WALK_COUNT = walks.size();

	Overlays overlays(WALK_COUNT, REF_COUNT);

	for (pt::u32 w_idx{}; w_idx < WALK_COUNT; ++w_idx) {

		std::set<pt::u32> walk_ref_idxs =
			find_refs_in_walk(g, walks.at(w_idx));

		for (auto &ref_idx : walk_ref_idxs) {
			const lq::ref_walk *ref_w =
				g.get_ref_vec(ref_idx)->walk;

			const pgt::walk_t &w = walks[w_idx];

			pt::u32 s_v_id = w.front().v_id;

			const std::vector<std::vector<pt::idx_t>> &s_vtx_refs =
				g.get_vertex_refs(s_v_id);

			const std::vector<pt::idx_t> &ref_idx_starts =
				s_vtx_refs[ref_idx];

			std::vector<overlay_prefix_sum_t> p_sums =
				overlay_walk(ref_w, w, ref_idx_starts);

			overlays.add_overlay(w_idx, ref_idx, std::move(p_sums));
		}
	}

	find_sub_inv(overlays, REF_COUNT, WALK_COUNT);

	return overlays;
}

std::set<pt::u32> find_refs_in_slice(const bd::VG &g, const pgt::walk_t &walk,
				     const pt::slice_t &sl,
				     pvr::var_type_e var_type)
{
	auto [start, len] = sl;
	auto [s_v_id, s_o] = walk[start];
	auto [t_v_id, t_o] = len == 0 ? walk[start] : walk[len - 1];

	const std::vector<std::vector<pt::idx_t>> &s_vtx_refs =
		g.get_vertex_refs(s_v_id);
	const std::vector<std::vector<pt::idx_t>> &t_vtx_refs =
		g.get_vertex_refs(t_v_id);

	std::set<pt::id_t> graph_walk_refs;

	for (pt::u32 ref_idx{}; ref_idx < g.ref_count(); ref_idx++) {
		if (s_vtx_refs[ref_idx].empty() || t_vtx_refs[ref_idx].empty())
			continue;

		graph_walk_refs.insert(ref_idx);
	}

	return graph_walk_refs;
}

void pop_exp(const bd::VG &g, const std::vector<pgt::walk_t> &walks,
	     const Overlays &ov, const pvr::raw_variant &pv, pt::u32 wa_idx,
	     pt::u32 wb_idx, std::map<pt::u32, pga::itn_t> &ref_map,
	     std::map<pt::idx_t, std::set<pt::idx_t>> &walk_to_refs,
	     bool &is_tangled)
{
	const pt::u32 REF_COUNT = g.ref_count();

	const auto &[sl_a, sl_b, vt] = pv; // slice a, slice b, var type

	for (pt::u32 ref_idx{}; ref_idx < REF_COUNT; ref_idx++) {
		const std::vector<overlay_prefix_sum_t> &wa_pss =
			ov.get_ps(wa_idx, ref_idx);
		for (const auto &[ref_start_idx, ps_f, ps_r] : wa_pss) {
			auto [start, len] = sl_a;

			// add context for subs and insertions
			pt::u32 i;
			pt::u32 N;
			pt::u32 as_walk_start_idx;
			pt::u32 as_ref_start_idx;
			pt::u32 as_len;

			switch (vt) {
			case pvr::var_type_e::del:
				i = start - 1;
				N = start + len;
				as_walk_start_idx = start - 1;
				as_ref_start_idx = ref_start_idx;
				as_len = len + 1 + 1;
				break;
			case pvr::var_type_e::sub:
				i = start - 1;
				N = start + len;
				as_walk_start_idx = start - 1;
				as_ref_start_idx = ref_start_idx;
				as_len = len + 1 + 1;
				break;
			case pvr::var_type_e::ins:
				i = start;
				N = start + 1; // len is 0 for insertions
				as_walk_start_idx = start;
				as_ref_start_idx = ref_start_idx;
				as_len = 1 + 1;
				break;
			}

			// Check prefix sums in one step
			bool is_forward_zero =
				(ps_f[N] - ps_f[i]) + ps_f[N] == 0;
			bool is_reverse_zero =
				(ps_r[N] - ps_r[i]) + ps_r[N] == 0;

			if (!is_forward_zero && !is_reverse_zero)
				continue;

			auto as_or = is_forward_zero ? pgt::or_e::forward
						     : pgt::or_e::reverse;

			if (as_or == pgt::or_e::forward) {
				// adjust ref start idx for reverse
				as_ref_start_idx = ref_start_idx + i;
			}
			else {
				as_ref_start_idx = ref_start_idx - i;
			}

			// valid overlay
			pga::allele_slice_t at{&walks.at(wa_idx),
					       wa_idx,
					       as_walk_start_idx,
					       g.get_ref_vec(ref_idx)->walk,
					       ref_idx,
					       as_ref_start_idx,
					       as_len,
					       as_or,
					       vt};

			pga::itn_t &itn = ref_map[ref_idx];
			itn.append_at(std::move(at));

			walk_to_refs[wa_idx].insert(ref_idx);
			if (itn.at_count() > 1)
				is_tangled = true;
		}
	}

	for (pt::u32 ref_idx{}; ref_idx < REF_COUNT; ref_idx++) {
		const std::vector<overlay_prefix_sum_t> &wb_pss =
			ov.get_ps(wb_idx, ref_idx);
		for (const auto &[ref_start_idx, ps_f, ps_r] : wb_pss) {
			auto [start, len] = sl_b;
			auto vt_ = pvr::covariant(vt);

			pt::u32 i;
			pt::u32 N;
			pt::u32 as_walk_start_idx;
			pt::u32 as_ref_start_idx;
			pt::u32 as_len;

			switch (vt_) {
			case pvr::var_type_e::del:
				i = start - 1;
				N = start + len;
				as_walk_start_idx = start - 1;
				as_ref_start_idx = ref_start_idx;
				as_len = len + 1 + 1;
				break;
			case pvr::var_type_e::sub:
				i = start - 1;
				N = start + len;
				as_walk_start_idx = start - 1;
				as_ref_start_idx = ref_start_idx;
				as_len = len + 1 + 1;
				break;
			case pvr::var_type_e::ins:
				i = start;
				N = start + 1; // len is 0 for insertions
				as_walk_start_idx = start;
				as_ref_start_idx = ref_start_idx;
				as_len = 1 + 1;
				break;
			}

			// Check prefix sums in one step
			bool is_forward_zero =
				(ps_f[N] - ps_f[i]) + ps_f[N] == 0;
			bool is_reverse_zero =
				(ps_r[N] - ps_r[i]) + ps_r[N] == 0;

			if (!is_forward_zero && !is_reverse_zero)
				continue;

			auto as_or = is_forward_zero ? pgt::or_e::forward
						     : pgt::or_e::reverse;

			if (as_or == pgt::or_e::forward) {
				// adjust ref start idx for reverse
				as_ref_start_idx = ref_start_idx + i;
			}
			else {
				as_ref_start_idx = ref_start_idx - i;
			}

			pga::allele_slice_t at{&walks.at(wb_idx),
					       wb_idx,
					       as_walk_start_idx,
					       g.get_ref_vec(ref_idx)->walk,
					       ref_idx,
					       as_ref_start_idx,
					       as_len,
					       as_or,
					       vt_};

			pga::itn_t &itn = ref_map[ref_idx];
			itn.append_at(std::move(at));

			walk_to_refs[wb_idx].insert(ref_idx);
			if (itn.at_count() > 1)
				is_tangled = true;
		}
	}
}

/**
 * [out] rov_exps: vector of expeditions, one per pairwise variant set
 */
std::pair<std::vector<pga::Exp>, std::vector<sub_inv>>
comp_overlays3(const bd::VG &g, const pvr::RoV &rov,
	       const std::set<pt::id_t> &to_call_ref_ids)
{
	std::vector<pga::Exp> rov_exps;

	const std::vector<pgt::walk_t> &walks = rov.get_walks();
	const std::vector<pvr::pairwise_variants> &pv = rov.get_irreducibles();

	// if (walks.size() > 500) {
	//	std::cerr << "Skipping " << rov.as_str() << " ---- "
	//		  << walks.size() << "----" << pv.size() << "\n";
	//	return {{}, {}};
	// }

#ifdef DEBUG
	if (pv.empty()) {
		ERR("No pairwise variants in RoV {}", rov.as_str());
		std::exit(EXIT_FAILURE);
	}
#endif

	Overlays ov = overlay_walks(g, walks);

	return {{}, {}};

	for (const pvr::pairwise_variants &p : pv) {

		auto [wa_idx, wb_idx, variants] = p;

		for (pt::u32 i{}; i < p.size(); i++) {
			pga::Exp e(&rov);
			std::map<pt::id_t, pga::itn_t> &ref_map =
				e.get_ref_itns_mut();
			std::map<pt::idx_t, std::set<pt::idx_t>> &walk_to_refs =
				e.get_walk_to_ref_idxs_mut();

			bool is_tangled = false;
			pop_exp(g, walks, ov, variants.at(i), wa_idx, wb_idx,
				ref_map, walk_to_refs, is_tangled);

			std::set<pt::id_t> exp_refs = e.get_ref_ids();

			// get set intersection with to_call_ref_ids
			std::set<pt::id_t> intersect_refs;
			std::set_intersection(
				exp_refs.begin(), exp_refs.end(),
				to_call_ref_ids.begin(), to_call_ref_ids.end(),
				std::inserter(intersect_refs,
					      intersect_refs.begin()));

			if (intersect_refs.empty())
				continue;

			e.set_tangled(is_tangled);

			rov_exps.emplace_back(std::move(e));
		}
	}

	return {rov_exps, {}};
}

std::pair<std::vector<pga::Exp>, std::vector<sub_inv>>
comp_itineraries3(const bd::VG &g, const pvr::RoV &rov,
		  const std::set<pt::id_t> &to_call_ref_ids)
{
	// auto can_be_non_planar = [](const pvst::vf_e fam) -> bool
	// {
	//	return fam == pvst::vf_e::flubble;
	// };

	if (!rov.can_be_non_planar())
		return overlay_tiny(g, rov, to_call_ref_ids);

	return generic::overlay_generic(g, rov, to_call_ref_ids);
}

} // namespace povu::overlay

namespace povu::overlay::shared
{
void update_exp(pt::u32 r_idx, pt::u32 w_idx, pga::allele_slice_t &&at,
		pga::Exp &e)
{
	std::map<pt::id_t, pga::itn_t> &ref_map = e.get_ref_itns_mut();
	pga::itn_t &itn = ref_map[r_idx];

	std::map<pt::idx_t, std::set<pt::idx_t>> &walk_to_refs =
		e.get_walk_to_ref_idxs_mut();

	itn.append_at_sorted(std::move(at));
	walk_to_refs[w_idx].insert(r_idx);
}
} // namespace povu::overlay::shared
