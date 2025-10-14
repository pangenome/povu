#include "povu/genomics/allele.hpp"

#include <csignal> // for raise, SIGINT

#include <cstdint>
#include <cstdio>
#include <cstdlib> // for exit, EXIT_FAILURE
#include <iostream>
#include <liteseq/refs.h> // for ref_walk, ref
#include <map>		  // for map
#include <utility>
#include <vector> // for vector

#include "povu/common/core.hpp"
#include "povu/common/log.hpp"
#include "povu/common/utils.hpp"
#include "povu/genomics/rov.hpp"
#include "povu/graph/types.hpp"

namespace povu::genomics::allele
{
namespace lq = liteseq;

constexpr pgr::var_type_e ins = pgr::var_type_e::ins;
constexpr pgr::var_type_e del = pgr::var_type_e::del;
constexpr pgr::var_type_e sub = pgr::var_type_e::sub;

const pt::u8 SLICE_A_IDX{0};
const pt::u8 SLICE_B_IDX{1};
using prefix_sum = std::vector<pt::u32>;
using rs_to_ps = std::pair<pt::u32, prefix_sum>; // ref start to prefix sum

class Overlays
{
	// i is the walk idx and j is the ref idx
	std::vector<std::vector<std::vector<rs_to_ps>>> prefix_sums_;
	// key is walk idx and value is idx in overlays vector
	// std::map<pt::u32, pt::u32> walk_to_overlay_;

public:
	Overlays(pt::u32 walk_count, pt::u32 ref_count)
	    : prefix_sums_(walk_count,
			   std::vector<std::vector<rs_to_ps>>(ref_count))
	{}

	[[nodiscard]]
	const std::vector<rs_to_ps> &get_ps(pt::u32 w_idx, pt::u32 r_idx) const
	{
		return prefix_sums_.at(w_idx).at(r_idx);
	}

	void add_overlay(pt::u32 w_idx, pt::u32 r_idx,
			 std::vector<rs_to_ps> &&p_sums)
	{
		prefix_sums_.at(w_idx).at(r_idx) = std::move(p_sums);
	}
};

struct overlay_t {
	pt::idx_t graph_w_start_idx;
	pt::idx_t ref_start_idx;
	pt::idx_t len;
	ptg::or_e slice_or;
	pgr::var_type_e vt;
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

std::vector<rs_to_ps> overlay_walk(const lq::ref_walk *ref_w,
				   const pgt::walk_t &graph_w,
				   const std::vector<pt::idx_t> &ref_idx_starts)
{

	const pt::u32 GRAPH_W_LEN = graph_w.size();
	std::vector<pt::u8> ov(GRAPH_W_LEN, 0);

	// std::string w_str = pgt::to_string(graph_w);
	// std::cerr << w_str << "\t";

	std::vector<rs_to_ps> prefix_sums;
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

		prefix_sums.emplace_back(ref_idx, comp_prefix_sum(ov));

		// reset ov for next ref start
		std::fill(ov.begin(), ov.end(), 0);
	}

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

/**
 * comp ps for each walk as maps to the ref
 */
Overlays comp_prefixes(const bd::VG &g, const std::vector<pgt::walk_t> &walks)
{
	const pt::u32 REF_COUNT = g.ref_count();
	const pt::u32 WALK_COUNT = walks.size();

	Overlays overlays(WALK_COUNT, REF_COUNT);

	for (pt::u32 w_idx{}; w_idx < WALK_COUNT; ++w_idx) {

		std::set<pt::u32> walk_ref_idxs =
			find_refs_in_walk(g, walks.at(w_idx));

		// std::cerr << pgt::to_string(walks.at(w_idx)) << "\n";
		// std::cerr << "refs: ";
		// povu::utils::print_with_comma(std::cerr, walk_ref_idxs, ',');
		// std::cerr << "\n";

		for (auto &ref_idx : walk_ref_idxs) {
			const lq::ref_walk *ref_w =
				g.get_ref_vec(ref_idx)->walk;

			const pgt::walk_t &w = walks[w_idx];

			pt::u32 s_v_id = w.front().v_id;

			const std::vector<std::vector<pt::idx_t>> &s_vtx_refs =
				g.get_vertex_refs(s_v_id);

			const std::vector<pt::idx_t> &ref_idx_starts =
				s_vtx_refs[ref_idx];

			std::vector<rs_to_ps> p_sums =
				overlay_walk(ref_w, w, ref_idx_starts);

			// if (dbg && ref_idx == 10) {
			//	std::cerr << "w idx " << w_idx << "\n";
			//	for (auto &[ref_start_idx, ps] : p_sums) {
			//		std::cerr << "ref idx " << ref_idx
			//			  << " ref start idx "
			//			  << ref_start_idx << " ps: ";
			//		povu::utils::print_with_comma(std::cerr,
			//					      ps, ',');
			//		std::cerr << "\n";
			//	}
			// }

			overlays.add_overlay(w_idx, ref_idx, std::move(p_sums));
		}
	}

	return overlays;
}

std::set<pt::u32> find_refs_in_slice(const bd::VG &g, const pgt::walk_t &walk,
				     const pt::slice_t &sl,
				     pgr::var_type_e var_type)
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
	     const Overlays &ov, const pgr::raw_variant &pv, pt::u32 wa_idx,
	     pt::u32 wb_idx, std::map<pt::u32, itn_t> &ref_map,
	     std::map<pt::idx_t, std::set<pt::idx_t>> &walk_to_refs,
	     bool &is_tangled)
{
	const pt::u32 REF_COUNT = g.ref_count();

	const auto &[sl_a, sl_b, vt] = pv; // slice a, slice b, var type

	for (pt::u32 ref_idx{}; ref_idx < REF_COUNT; ref_idx++) {
		const std::vector<rs_to_ps> &wa_pss =
			ov.get_ps(wa_idx, ref_idx);
		for (const auto &[ref_start_idx, ps] : wa_pss) {
			auto [start, len] = sl_a;

			// add context for subs and insertions
			pt::u32 i;
			pt::u32 N;
			pt::u32 as_walk_start_idx;
			pt::u32 as_ref_start_idx;
			pt::u32 as_len;

			switch (vt) {
			case pgr::var_type_e::del:
				i = start - 1;
				N = start + len;
				as_walk_start_idx = start - 1;
				as_ref_start_idx = ref_start_idx;
				as_len = len + 1 + 1;
				break;
			case pgr::var_type_e::sub:
				i = start - 1;
				N = start + len;
				as_walk_start_idx = start - 1;
				as_ref_start_idx = ref_start_idx;
				as_len = len + 1 + 1;
				break;
			case pgr::var_type_e::ins:
				i = start;
				N = start + 1; // len is 0 for insertions
				as_walk_start_idx = start;
				as_ref_start_idx = ref_start_idx;
				as_len = 1 + 1;
				break;
			}

			if ((ps[N] - ps[i]) == 0) {
				// valid overlay
				allele_slice_t at{&walks.at(wa_idx),
						  wa_idx,
						  as_walk_start_idx,
						  g.get_ref_vec(ref_idx)->walk,
						  ref_idx,
						  as_ref_start_idx,
						  as_len,
						  pgt::or_e::forward,
						  vt};

				itn_t &itn = ref_map[ref_idx];
				itn.append_at(std::move(at));
				walk_to_refs[wa_idx].insert(ref_idx);
				if (itn.at_count() > 1)
					is_tangled = true;
			}
		}
	}

	for (pt::u32 ref_idx{}; ref_idx < REF_COUNT; ref_idx++) {
		const std::vector<rs_to_ps> &wb_pss =
			ov.get_ps(wb_idx, ref_idx);
		for (const auto &[ref_start_idx, ps] : wb_pss) {
			auto [start, len] = sl_b;
			auto vt_ = pgr::covariant(vt);

			pt::u32 i;
			pt::u32 N;
			pt::u32 as_walk_start_idx;
			pt::u32 as_ref_start_idx;
			pt::u32 as_len;

			switch (vt_) {
			case pgr::var_type_e::del:
				i = start - 1;
				N = start + len;
				as_walk_start_idx = start - 1;
				as_ref_start_idx = ref_start_idx;
				as_len = len + 1 + 1;
				break;
			case pgr::var_type_e::sub:
				i = start - 1;
				N = start + len;
				as_walk_start_idx = start - 1;
				as_ref_start_idx = ref_start_idx;
				as_len = len + 1 + 1;
				break;
			case pgr::var_type_e::ins:
				i = start;
				N = start + 1; // len is 0 for insertions
				as_walk_start_idx = start;
				as_ref_start_idx = ref_start_idx;
				as_len = 1 + 1;
				break;
			}

			if ((ps[N] - ps[i]) == 0) { // valid overlay
				allele_slice_t at{&walks.at(wb_idx),
						  wb_idx,
						  as_walk_start_idx,
						  g.get_ref_vec(ref_idx)->walk,
						  ref_idx,
						  as_ref_start_idx,
						  as_len,
						  pgt::or_e::forward,
						  vt_};

				itn_t &itn = ref_map[ref_idx];
				itn.append_at(std::move(at));
				walk_to_refs[wb_idx].insert(ref_idx);
				if (itn.at_count() > 1)
					is_tangled = true;
			}
		}
	}
}

/**
 * [out] rov_exps: vector of expeditions, one per pairwise variant set
 */
std::vector<Exp> comp_overlays3(const bd::VG &g, const pgr::RoV &rov)
{
	std::vector<Exp> rov_exps;

	const std::vector<pgt::walk_t> &walks = rov.get_walks();
	const std::vector<pgr::pairwise_variants> &pv = rov.get_irreducibles();

#ifdef DEBUG
	if (pv.empty()) {
		ERR("No pairwise variants in RoV {}", rov.as_str());
		std::exit(EXIT_FAILURE);
	}
#endif

	Overlays ov = comp_prefixes(g, walks);

	for (const rov::pairwise_variants &p : pv) {

		auto [wa_idx, wb_idx, variants] = p;

		for (pt::u32 i{}; i < p.size(); i++) {
			Exp e(&rov);
			std::map<pt::id_t, itn_t> &ref_map =
				e.get_ref_itns_mut();
			std::map<pt::idx_t, std::set<pt::idx_t>> &walk_to_refs =
				e.get_walk_to_ref_idxs_mut();

			bool is_tangled = false;
			pop_exp(g, walks, ov, variants.at(i), wa_idx, wb_idx,
				ref_map, walk_to_refs, is_tangled);

			e.set_tangled(is_tangled);

			rov_exps.emplace_back(std::move(e));
		}
	}

	return rov_exps;
}

std::vector<Exp> comp_itineraries3(const bd::VG &g, const pgr::RoV &rov)
{

	return comp_overlays3(g, rov);
}

} // namespace povu::genomics::allele
