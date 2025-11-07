#include "povu/genomics/allele.hpp"

#include <csignal> // for raise, SIGINT

#include <cstdint>
#include <cstdio>
#include <cstdlib> // for exit, EXIT_FAILURE
#include <iostream>
#include <liteseq/refs.h> // for ref_walk, ref
#include <map>		  // for map
#include <vector>	  // for vector

#include "povu/common/core.hpp"
#include "povu/common/log.hpp"
#include "povu/common/utils.hpp"
#include "povu/graph/types.hpp"
#include "povu/variation/rov.hpp"

namespace povu::genomics::allele
{
namespace lq = liteseq;

constexpr pvr::var_type_e ins = pvr::var_type_e::ins;
constexpr pvr::var_type_e del = pvr::var_type_e::del;
constexpr pvr::var_type_e sub = pvr::var_type_e::sub;

bool operator!=(const allele_slice_t &lhs, const allele_slice_t &rhs)
{
	return !(lhs == rhs);
}

bool operator==(const allele_slice_t &lhs, const allele_slice_t &rhs)
{
	if (lhs.len != rhs.len)
		return false;

	pt::u32 N = lhs.len;
	for (pt::u32 i{}; i < N; i++) {
		auto [lhs_v_id, lhs_o] = (*lhs.walk)[lhs.walk_start_idx + i];
		auto [rhs_v_id, rhs_o] = (*rhs.walk)[rhs.walk_start_idx + i];

		if (lhs_v_id != rhs_v_id || lhs_o != rhs_o)
			return false;
	}

	return true;
}

bool is_contained(const std::vector<pt::slice_t> &ref_slices,
		  pt::slice_t ref_slice)
{
	for (const pt::slice_t &s : ref_slices) {
		// check if is prefix
		if (s.start == ref_slice.start && s.len > ref_slice.len) {
			return true;
		}

		// check if is suffix
		if (s.start < ref_slice.start &&
		    s.start + s.len == ref_slice.start + ref_slice.len) {
			return true;
		}

		// check if is contained
		if (s.start > ref_slice.start &&
		    s.start + s.len < ref_slice.start + ref_slice.len) {
			return true;
		}
	}

	return false;
}

pt::idx_t is_valid_rev(const lq::ref_walk *ref_w, pt::idx_t ref_w_start_idx,
		       const pgt::walk_t &w, pt::idx_t w_start_idx,
		       pt::idx_t len)
{
	pt::idx_t valid_len{0};

	for (pt::idx_t i{}; i < len; i++) {
		pt::idx_t ref_w_idx = ref_w_start_idx - i;
		if (ref_w_idx > ref_w->step_count)
			return valid_len;

		pt::idx_t ref_v_id = ref_w->v_ids[ref_w_idx];
		pgt::or_e ref_o =
			ref_w->strands[ref_w_idx] == lq::strand::STRAND_FWD
				? pgt::or_e::forward
				: pgt::or_e::reverse;

		pt::idx_t graph_w_idx = w_start_idx + i;
		auto [w_v_id, w_o] = w[graph_w_idx];

		if (ref_v_id != w_v_id || ref_o != ptg::flip(w_o)) {
			break;
		}

		valid_len++;
	}

	return valid_len;
}

pt::idx_t is_valid(const lq::ref_walk *ref_w, pt::idx_t ref_w_start_idx,
		   const pgt::walk_t &w, pt::idx_t w_start_idx, pt::idx_t len)
{
	pt::idx_t valid_len{0};

	for (pt::idx_t i{}; i < len; i++) {
		pt::idx_t ref_w_idx = ref_w_start_idx + i;
		pt::idx_t graph_w_idx = w_start_idx + i;

		pt::idx_t ref_v_id = ref_w->v_ids[ref_w_idx];
		pgt::or_e ref_o =
			ref_w->strands[ref_w_idx] == lq::strand::STRAND_FWD
				? pgt::or_e::forward
				: pgt::or_e::reverse;

		//  auto [ref_v_id, ref_o, _] = ref_w[ref_w_idx];
		auto [w_v_id, w_o] = w[graph_w_idx];

		if (ref_v_id != w_v_id || ref_o != w_o)
			break;

		valid_len++;
	}

	return valid_len;
}

void run_conv(pt::idx_t ref_idx, const lq::ref_walk *ref_w,
	      std::map<pt::idx_t, std::set<pt::idx_t>> &walk_to_refs,
	      pt::idx_t w_idx,
	      const std::pair<const pgt::walk_t &, pt::slice_t> &walk_slice,
	      std::vector<pt::slice_t> &walk_slices, itn_t &ref_itns,
	      const std::vector<pt::idx_t> &vtx_ref_idxs, bool &is_tangled)
{
	const pgt::walk_t &w = walk_slice.first;
	pt::idx_t walk_step_idx = walk_slice.second.start;
	pt::idx_t slice_len = walk_slice.second.len;

	for (pt::idx_t ref_step_idx : vtx_ref_idxs) {
		// Assumption:
		// in a valid GFA file,
		// if a ref starts within a walk then it has to
		// start at the beginning of the walk
		if (walk_step_idx > 0 && ref_step_idx != 0)
			continue;

		std::pair<const lq::ref_walk *, pt::slice_t> ref_slice = {
			ref_w, {ref_step_idx, slice_len}};

		ptg::or_e ref_or = pgt::or_e::forward;
		pt::idx_t valid_len{};

		valid_len = is_valid(ref_w, ref_step_idx, w, walk_step_idx,
				     slice_len);

		if (valid_len >= 2) {
			ref_or = pgt::or_e::forward;
		}
		else { // valid_len < 2
			valid_len = is_valid_rev(ref_w, ref_step_idx, w,
						 walk_step_idx, slice_len);
			ref_or = pgt::or_e::reverse;
		}

		if (valid_len < 2)
			continue;

		if (!is_contained(walk_slices, {walk_step_idx, valid_len})) {
			walk_slices.push_back({walk_step_idx, valid_len});
			ref_itns.append_at({&w, w_idx, walk_step_idx, ref_w,
					    ref_idx, ref_step_idx, valid_len,
					    ref_or});
			walk_to_refs[w_idx].insert(ref_idx);

			if (ref_itns.at_count() > 1)
				is_tangled = true;
		}
	}
}

bool comp_overlays(const bd::VG &g, const pgt::walk_t &w, pt::idx_t w_idx,
		   std::map<pt::id_t, itn_t> &ref_map,
		   std::map<pt::idx_t, std::set<pt::idx_t>> &walk_to_refs)
{
	bool is_tangled{false};
	const pt::idx_t W_LEN = w.size();

	for (pt::idx_t ref_idx{}; ref_idx < g.get_ref_count(); ref_idx++) {
		const lq::ref_walk *ref_w = g.get_ref_vec(ref_idx)->walk;
		itn_t ref_itns;
		std::vector<pt::slice_t> walk_slices;

		for (pt::idx_t w_step_idx{}; w_step_idx < W_LEN; ++w_step_idx) {
			pt::idx_t slice_len = W_LEN - w_step_idx;
			const pgt::step_t &step = w[w_step_idx];
			auto [v_id, o] = step;
			pt::idx_t v_idx = g.v_id_to_idx(v_id);
			const std::vector<pt::idx_t> &vtx_ref_idxs =
				g.get_vertex_ref_idxs(v_idx, ref_idx);

			std::pair<const pgt::walk_t &, pt::slice_t> walk_slice =
				{w, {w_step_idx, slice_len}};

			run_conv(ref_idx, ref_w, walk_to_refs, w_idx,
				 walk_slice, walk_slices, ref_itns,
				 vtx_ref_idxs, is_tangled);
		}

		if (ref_itns.at_count() > 0)
			ref_map.emplace(ref_idx, std::move(ref_itns));
	}

	return is_tangled;
}

bool overlay_leftwards(const lq::ref_walk *ref_w, const pgt::walk_t &graph_w,
		       const std::vector<pt::idx_t> &vtx_ref_idxs,
		       pt::u32 ref_w_start_idx, pt::u32 graph_w_start_idx,
		       pt::u32 len)
{

	pt::u32 valid_len{0};
	for (pt::u32 i{}; i < len; i++) {
		pt::idx_t ref_w_idx = ref_w_start_idx + i;
		pt::idx_t graph_w_idx = graph_w_start_idx + i;

		pt::idx_t ref_v_id = ref_w->v_ids[ref_w_idx];
		pgt::or_e ref_o =
			ref_w->strands[ref_w_idx] == lq::strand::STRAND_FWD
				? pgt::or_e::forward
				: pgt::or_e::reverse;

		//  auto [ref_v_id, ref_o, _] = ref_w[ref_w_idx];
		auto [w_v_id, w_o] = graph_w[graph_w_idx];

		if (ref_v_id != w_v_id || ref_o != w_o)
			return false;

		valid_len++;
	}
	return true;
}

struct overlay_t {
	pt::idx_t graph_w_start_idx;
	pt::idx_t ref_start_idx;
	pt::idx_t len;
	ptg::or_e slice_or;
	pvr::var_type_e vt;
};

const pt::u8 SLICE_A_IDX{0};
const pt::u8 SLICE_B_IDX{1};

/**
 * [out] walk_to_refs: map of walk idx to ref idxs that take the walk
 */
std::map<pt::u32, std::vector<overlay_t>>
overlay(const bd::VG &g, const pgt::walk_t &graph_w,
	const std::set<pt::u32> &graph_walk_refs,
	const std::vector<pvr::raw_variant> &variants, pt::u8 sl_idx, bool dbg)
{
	std::map<pt::u32, std::vector<overlay_t>> ref_to_overlays;

	std::string w_str = pgt::to_string(graph_w);

	const pt::u32 GRAPH_W_LEN = graph_w.size();
	const pt::u32 REF_COUNT = g.get_ref_count();

	for (pt::u32 ref_idx : graph_walk_refs) {
		const lq::ref_walk *ref_w = g.get_ref_vec(ref_idx)->walk;

		for (const pvr::raw_variant &v : variants) {

			auto [sl_a, sl_b, vt_] = v;
			auto [start, len] = sl_idx == SLICE_A_IDX ? sl_a : sl_b;
			pvr::var_type_e vt = sl_idx == SLICE_A_IDX
						     ? vt_
						     : pvr::covariant(vt_);
			pvr::var_type_e overlay_vt = vt;

			pt::u32 i{start};
			pt::u32 N = start + len;
			if (sl_idx == SLICE_A_IDX) {
				if (vt == sub || vt == ins) {
					N++;
					i--;
				}
				else if (vt == del) {
					N++;
				}
			}
			else { // SLICE_B_IDX
				if (vt == sub) {
					N++;
					i--;
				}
				else if (vt == ins) {
					N++;
				}
				else if (vt == del) {
					N++;
					i--;
				}
			}

			pt::idx_t slice_len = N - i;

			if (i > GRAPH_W_LEN)
				continue;

			auto [v_id, o] = graph_w.at(i);
			pt::idx_t v_idx = g.v_id_to_idx(v_id);

			const std::vector<pt::idx_t> &vtx_ref_idxs =
				g.get_vertex_ref_idxs(v_idx, ref_idx);

			if (vtx_ref_idxs.empty()) {
				continue;
			}

			for (pt::u32 ref_w_start_idx : vtx_ref_idxs) {
				bool a = overlay_leftwards(
					ref_w, graph_w, vtx_ref_idxs,
					ref_w_start_idx, i, slice_len);

				if (a) {
					overlay_t o{
						i, ref_w_start_idx, slice_len,
						pgt::or_e::forward, overlay_vt};

					ref_to_overlays[ref_idx].push_back(o);
				}
			}
		}
	}

	return ref_to_overlays;
}

/**
 * refs are in a walk if
 */
std::set<pt::id_t> refs_in_walk(const bd::VG &g, const pgt::walk_t &walk)
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

/**
 * [out] rov_exps: vector of expeditions, one per pairwise variant set
 */
void comp_overlays2(const bd::VG &g, const std::vector<pgt::walk_t> &walks,
		    const std::vector<pvr::pairwise_variants> &pv,
		    const pvr::RoV *rov, std::vector<Exp> &rov_exps)
{
	// bool dbg = rov->as_str() == ">1546>1551" ? true : false;

	bool dbg = walks.size() > 5 ? true : false;

	if (dbg) {
	}

	auto foo = [&](const std::map<pt::u32, std::vector<overlay_t>> &x,
		       pt::u32 w_idx, std::map<pt::id_t, itn_t> &ref_map,
		       std::map<pt::idx_t, std::set<pt::idx_t>> &walk_to_refs)
		-> bool
	{
		bool is_tangled{false};
		for (const auto &[ref_idx, overlays] : x) {
			if (overlays.size() > 1)
				is_tangled = true;

			walk_to_refs[w_idx].insert(ref_idx);
			itn_t &itn = ref_map[ref_idx];
			for (const auto &o : overlays) {
				itn.append_at(allele_slice_t{
					&walks.at(w_idx), w_idx,
					o.graph_w_start_idx,
					g.get_ref_vec(ref_idx)->walk, ref_idx,
					o.ref_start_idx, o.len, o.slice_or,
					o.vt});
			}
		}

		return is_tangled;
	};

	if (dbg)
		volatile int z = 0;

	if (dbg) {
		std::cerr << "Computing overlays for RoV " << rov->as_str()
			  << " with " << walks.size() << " walks and "
			  << pv.size() << " pairwise variants\n";

		std::cerr << "Walks for RoV " << rov->as_str() << ":\n";
		for (int i = 0; i < walks.size(); ++i) {
			std::cerr << i << ": " << pgt::to_string(walks[i])
				  << "\n";
		}
	}

	pt::Time::time_point start;
	pt::Time::time_point end;

	start = pt::Time::now();

	std::map<pt::u32, std::set<pt::id_t>> refs_in_walks;
	for (pt::u32 w_idx{}; w_idx < walks.size(); ++w_idx) {
		refs_in_walks[w_idx] = refs_in_walk(g, walks[w_idx]);
	}

	end = pt::Time::now();
	auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(
		end - start);

	if (dbg) {
		std::cerr << __func__
			  << " Elapsed time (1): " << elapsed.count() << " ns"
			  << "\n";
	}
	start = pt::Time::now();

	for (const pvr::pairwise_variants &p : pv) {
		Exp e(rov);
		std::map<pt::id_t, itn_t> &ref_map = e.get_ref_itns_mut();
		std::map<pt::idx_t, std::set<pt::idx_t>> &walk_to_refs =
			e.get_walk_to_ref_idxs_mut();

		auto [w1_idx, w2_idx, variants] = p;

		if (dbg) {
			std::cerr << "(" << w1_idx << "," << w2_idx << ")\n";
			std::cerr << "Variants:\n";
			for (const pvr::raw_variant &rv : variants) {
				std::cerr << rv << "\n";
			}
		}

		auto x = overlay(g, walks.at(w1_idx), refs_in_walks.at(w1_idx),
				 variants, SLICE_A_IDX, dbg);

		bool is_w1_tangled = foo(x, w1_idx, ref_map, walk_to_refs);

		auto y = overlay(g, walks.at(w2_idx), refs_in_walks.at(w1_idx),
				 variants, SLICE_B_IDX, dbg);

		bool is_w2_tangled = foo(y, w2_idx, ref_map, walk_to_refs);

		if (is_w1_tangled || is_w2_tangled)
			e.set_tangled(true);

		rov_exps.emplace_back(std::move(e));
	}

	// Calculate the elapsed time
	end = pt::Time::now();
	elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(end -
								       start);

	// Output the elapsed time in nanoseconds
	if (dbg)
		std::cerr << __func__
			  << " Elapsed time (2): " << elapsed.count() << " ns"
			  << "\n";

	if (">2597>2621" == rov->as_str()) {
		std::exit(1);
	}

	return;
}

std::vector<Exp> comp_itineraries2(const bd::VG &g, const pvr::RoV &rov)
{

	const std::vector<pgt::walk_t> &walks = rov.get_walks();
	const std::vector<pvr::pairwise_variants> &pv = rov.get_irreducibles();
	const pvr::RoV *rov_ = &rov;

	if (pv.empty()) {
		ERR("No pairwise variants in RoV {}", rov.as_str());
		std::exit(EXIT_FAILURE);
	}

	std::vector<Exp> rov_exps;
	comp_overlays2(g, walks, pv, rov_, rov_exps);
	return rov_exps;
}

void comp_itineraries(const bd::VG &g, Exp &exp)
{
	if (exp.get_rov() == nullptr) {
		ERR("RoV pointer is null");
		std::exit(EXIT_FAILURE);
	}

	const std::vector<pgt::walk_t> &walks = exp.get_rov()->get_walks();
	std::map<pt::id_t, itn_t> &ref_map = exp.get_ref_itns_mut();
	std::map<pt::idx_t, std::set<pt::idx_t>> &walk_to_refs =
		exp.get_walk_to_ref_idxs_mut();

	const std::vector<pvr::pairwise_variants> &pv =
		exp.get_rov()->get_irreducibles();

	for (pt::idx_t w_idx = 0; w_idx < walks.size(); ++w_idx) {

		bool is_tangled = comp_overlays(g, walks.at(w_idx), w_idx,
						ref_map, walk_to_refs);
		if (is_tangled)
			exp.set_tangled(true);

#ifdef DEBUG
		// ensure walk pointers are valid
		auto refs = exp.get_ref_idxs_for_walk(w_idx);
		for (auto r : refs) {
			const itn_t &r_itn = exp.get_itn(r);
			for (pt::idx_t at_idx = 0; at_idx < r_itn.at_count();
			     ++at_idx) {
				const allele_slice_t &as = r_itn.get_at(at_idx);
				auto w_idx = as.walk_idx;
				const pgt::walk_t &w =
					exp.get_rov()->get_walks().at(w_idx);
				if (as.walk != &w) {
					ERR("Inconsistent walk pointer in "
					    "allele_slice, {}",
					    exp.id());
					std::exit(EXIT_FAILURE);
				}

				pt::idx_t step_idx = as.walk_start_idx;
				pt::idx_t end = as.walk_start_idx + as.len;
				for (step_idx; step_idx < end; ++step_idx)
					auto _ = as.get_step(step_idx);
			}
		}
#endif
	}

	return;
}

} // namespace povu::genomics::allele
