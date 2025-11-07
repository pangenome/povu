#include "povu/overlay/sne.hpp"

// std includes
#include <algorithm> // for min, max
#include <optional>
#include <ostream>

#include <set>
#include <utility>
#include <vector>

// deps
#include <liteseq/refs.h> // for ref_walk, ref

// povu includes
#include "povu/common/constants.hpp"
#include "povu/common/core.hpp"
#include "povu/genomics/allele.hpp"
#include "povu/graph/pvst.hpp"
#include "povu/graph/types.hpp" // pgt
#include "povu/overlay/shared.hpp"
#include "povu/variation/rov.hpp"

namespace povu::overlay::sne
{
namespace lq = liteseq;
namespace pgt = povu::types::graph;

using pga::ref_needle;
using pga::sub_inv;

constexpr pgt::or_e fo = pgt::or_e::forward;
constexpr pgt::or_e ro = pgt::or_e::reverse;

/**
 * @brief match ref walks at index i and j in ref walk 1 and ref walk 2
 *
 * i index in ref walk 1
 * j index in ref walk 2
 * dist accumulative length of matched region
 */
pt::status_t do_match(const lq::ref_walk *ref_w1, const lq::ref_walk *ref_w2,
		      pt::u32 &i, pt::u32 &j, pt::u32 &dist)
{
	if (i == pc::INVALID_IDX || j == pc::INVALID_IDX)
		return -1;

	pt::idx_t v_id1 = ref_w1->v_ids[i];
	pgt::or_e o1 = ref_w1->strands[i] == lq::strand::STRAND_FWD
			       ? pgt::or_e::forward
			       : pgt::or_e::reverse;

	pt::idx_t v_id2 = ref_w2->v_ids[j];
	pgt::or_e o2 = ref_w2->strands[j] == lq::strand::STRAND_FWD
			       ? pgt::or_e::forward
			       : pgt::or_e::reverse;

	// compare v_id1,o1 with v_id2,o2
	if (v_id1 != v_id2 || o1 != pgt::flip(o2))
		return 1;

	dist++;

	return 0;
}

/**
 * are adjacent extensions
 */
bool is_adj_or_overlap(const extension &curr_ext, const extension &next_ext)
{
	bool x = curr_ext.f_ref_start < next_ext.f_ref_start &&
		 (curr_ext.f_ref_start + curr_ext.len) >= next_ext.f_ref_start;

	bool y = curr_ext.r_ref_start > next_ext.r_ref_start &&
		 (curr_ext.r_ref_start - curr_ext.len) <= next_ext.r_ref_start;

	return x && y;
}

bool is_contained(const extension &curr_ext, const extension &next_ext)
{
	bool x = curr_ext.f_ref_start < next_ext.f_ref_start &&
		 curr_ext.f_ref_start + curr_ext.len >=
			 next_ext.f_ref_start + next_ext.len;

	bool y = curr_ext.r_ref_start > next_ext.r_ref_start &&
		 curr_ext.r_ref_start - curr_ext.len <=
			 next_ext.r_ref_start - next_ext.len;

	return x && y;
}

std::set<pt::u32> find_containments(const std::vector<extension> &extensions)
{
	std::set<pt::u32> to_remove;

	pt::u32 N = extensions.size();
	for (pt::u32 i = 0; i < N; i++) {
		for (pt::u32 j = 0; j < N; j++) {
			if (i == j)
				continue;

			if (is_contained(extensions[i], extensions[j]))
				to_remove.insert(j);
		}
	}

	return to_remove;
}

std::optional<pt::op_t<pt::u32>>
is_overlap(const std::vector<extension> &extensions, pt::u32 curr_idx,
	   pt::u32 next_idx)
{
	bool a = is_adj_or_overlap(extensions[curr_idx], extensions[next_idx]);
	bool b = is_adj_or_overlap(extensions[next_idx], extensions[curr_idx]);

	if (a)
		return pt::op_t<pt::u32>{curr_idx, next_idx};

	if (b)
		return pt::op_t<pt::u32>{next_idx, curr_idx};

	return std::nullopt;
}

std::vector<extension> filter_overlaps(const std::vector<extension> &extensions)
{
	std::set<pt::u32> contained = find_containments(extensions);

	std::vector<extension> filtered;
	for (pt::u32 i{}; i < extensions.size(); ++i) {
		if (pv_cmp::contains(contained, i))
			continue;

		filtered.push_back(extensions[i]);
	}

	return filtered;
}

extension stitch(const std::vector<extension> &extensions,
		 std::vector<pt::u32> &stitch_path)
{
	extension stitched;
	stitched.f_r_idx = extensions[0].f_r_idx;
	stitched.r_r_idx = extensions[0].r_r_idx;
	stitched.f_ref_start = extensions[stitch_path.front()].f_ref_start;
	stitched.r_ref_start = extensions[stitch_path.front()].r_ref_start;

	stitched.len = 0;

	stitched.len = extensions[stitch_path.back()].f_ref_start +
		       extensions[stitch_path.back()].len -
		       stitched.f_ref_start;

	// add context
	stitched.f_ref_start--;
	stitched.f_ref_start++;
	stitched.len++;

	return stitched;
}

// merge overlapping extensions
std::vector<extension>
merge_extensions(const std::vector<extension> &extensions, const bd::VG &g)
{
	// std::vector<extension> filtered = extensions;
	std::vector<extension> filtered = filter_overlaps(extensions);

	StitchGraph sg;
	std::set<pt::op_t<pt::u32>> overlaps;
	for (pt::u32 i{}; i < filtered.size() - 1; ++i) {
		if (auto opt_ov = is_overlap(filtered, i, i + 1)) {
			overlaps.insert(*opt_ov);
			auto [a, b] = *opt_ov;
			sg.add_edge(a, b);
		}
	}

	sg.compute_bfs_forest();

	std::vector<extension> merged;
	const std::vector<Tree> &stitch_forest = sg.get_forest();

	for (const auto &t : stitch_forest) {
		std::set<pt::u32> leaves = t.get_leaves();
		for (auto l : leaves) {
			std::vector<pt::u32> sp = t.get_stitch_path(l);
			extension s_ext = stitch(filtered, sp);
			merged.push_back(s_ext);
		}
	}

	for (auto i : sg.get_isolated())
		merged.push_back(filtered[i]);

	return merged;
}

/**
 * Adds 1 before division for ceiling behavior
 */
pt::u32 half_diff_ceil(pt::u32 a, pt::u32 b)
{
	pt::u32 diff = std::max(a, b) - std::min(a, b);
	return (diff + 1) / 2;
}

pga::Exp sne_exp()
{
	auto *v = new pvst::SnE();
	v->set_height(0);
	auto *r = new pvr::RoV(v);
	auto e = pga::Exp(r);

	return e;
}

void comp_expeditions(const bd::VG &g, const std::vector<extension> &extensions,
		      std::vector<pga::Exp> &exps)
{
	for (pt::u32 i{}; i < extensions.size(); ++i) {
		const auto &ext = extensions[i];
		auto [f_r_idx, r_r_idx, f_ref_start, r_ref_start, len] = ext;

		pga::Exp e = sne_exp();

		//  create allele slice for fwd
		pga::allele_slice_t at_f{nullptr,
					 0,
					 0,
					 g.get_ref_vec(f_r_idx)->walk,
					 f_r_idx,
					 f_ref_start,
					 len,
					 pgt::or_e::forward,
					 pvr::var_type_e::sub};

		shared::update_exp(f_r_idx, 0, std::move(at_f), e);

		// create allele slice for rev
		pga::allele_slice_t at_r{nullptr,
					 1,
					 0,
					 g.get_ref_vec(r_r_idx)->walk,
					 r_r_idx,
					 r_ref_start,
					 len,
					 pgt::or_e::reverse,
					 pvr::var_type_e::sub};

		shared::update_exp(r_r_idx, 1, std::move(at_r), e);

		exps.push_back(e);
	}
}

std::vector<sub_inv> select(const std::vector<pt::u32> &selected_idxs,
			    const std::vector<sub_inv> &anchors)
{
	std::vector<sub_inv> selected;
	// convert selected idxs to sets
	std::set<pt::u32> idx_set;
	for (auto idx : selected_idxs)
		idx_set.insert(idx);

	for (pt::u32 i{}; i < anchors.size(); ++i)
		if (pv_cmp::contains(idx_set, i))
			selected.push_back(anchors[i]);

	return selected;
}

std::set<pt::u32> find_alt_ref_ids(const std::vector<pin_cushion> &pcushions,
				   const std::set<pt::u32> &to_call_ref_ids)
{
	std::set<pt::u32> alt_ref_ids;
	for (const pin_cushion &pc : pcushions)
		for (auto r_idx : pc.get_refs())
			alt_ref_ids.insert(r_idx);

	for (pt::u32 ref_r_id : to_call_ref_ids)
		alt_ref_ids.erase(ref_r_id);

	return alt_ref_ids;
}

chain_t gen_colinear_chain(pt::u32 ref_r_id, pt::u32 alt_r_id,
			   const std::vector<pin_cushion> &pcushions)
{
	chain_t chain;
	for (const pin_cushion &pc : pcushions) {
		auto opt_ref_dir = pc.get_ref_direction(ref_r_id);
		auto opt_alt_dir = pc.get_ref_direction(alt_r_id);

		if (!opt_alt_dir || !opt_ref_dir)
			continue;

		if (*opt_ref_dir == *opt_alt_dir)
			continue;

		pin ref_pin = pc.get_ref_pin(ref_r_id).value();
		pin alt_pin = pc.get_ref_pin(alt_r_id).value();

		chain.links.emplace_back(chain_link::from_pin(ref_pin),
					 chain_link::from_pin(alt_pin));

		chain.len_++;
	}

	return chain;
}

void set_chain_link_limits(chain_t &c)
{
	for (pt::u32 i{}; i < c.len(); ++i) {
		pt::u32 limit_left = pc::INVALID_IDX;
		pt::u32 limit_right = pc::INVALID_IDX;

		if (i == 0 && i + 1 < c.links.size()) {
			auto idx_next = c.links[i + 1].ref_link.idx;
			auto idx_curr = c.links[i].ref_link.idx;
			limit_right = half_diff_ceil(idx_curr, idx_next);
		}
		else if (i == c.links.size() - 1 && i > 0) {
			auto idx_prev = c.links[i - 1].ref_link.idx;
			auto idx_curr = c.links[i].ref_link.idx;
			limit_left = half_diff_ceil(idx_curr, idx_prev);
		}
		else if (i > 0 && i + 1 < c.links.size()) {
			auto idx_prev = c.links[i - 1].ref_link.idx;
			auto idx_next = c.links[i + 1].ref_link.idx;
			auto idx_curr = c.links[i].ref_link.idx;
			limit_left = half_diff_ceil(idx_curr, idx_prev);
			limit_right = half_diff_ceil(idx_curr, idx_next);
		}

		c.links[i].ref_link.limit_left = limit_left;
		c.links[i].ref_link.limit_right = limit_right;

		c.links[i].alt_link.limit_left = limit_left;
		c.links[i].alt_link.limit_right = limit_right;
	}
}

std::optional<extension> extend_link(const bd::VG &g, const link_pair &lp,
				     pt::u32 lidx, pt::u32 ref_r_id,
				     pt::u32 alt_r_id)
{
	auto [ref_link, alt_link] = lp;

	auto [fwd_link, rev_link] =
		(ref_link.dir == hap_traversal::left)
			? std::make_pair(ref_link, alt_link)
			: std::make_pair(alt_link, ref_link);

	auto [f_r_id, r_r_id] = (ref_link.dir == hap_traversal::left)
					? std::make_pair(ref_r_id, alt_r_id)
					: std::make_pair(alt_r_id, ref_r_id);

	const lq::ref_walk *ref_w1 = g.get_ref_vec(ref_r_id)->walk;
	const lq::ref_walk *ref_w2 = g.get_ref_vec(alt_r_id)->walk;

	const lq::ref_walk *f_ref_w = g.get_ref_vec(f_r_id)->walk;
	const lq::ref_walk *r_ref_w = g.get_ref_vec(r_r_id)->walk;

	auto [f_ref_start, f_limit_left, f_limit_right, _] = fwd_link;
	auto [r_ref_start, r_limit_left, r_limit_right, __] = rev_link;

	pt::u32 i{f_ref_start};
	pt::u32 j{r_ref_start};
	pt::u32 dist_a{}; // distance extended

	// i is ref; j is alt
	pt::u32 S_i = (f_limit_left == pc::INVALID_IDX)
			      ? f_ref_start
			      : f_ref_start - f_limit_left - 1;
	pt::u32 S_j = (f_limit_right == pc::INVALID_IDX)
			      ? r_ref_start
			      : r_ref_start - r_limit_right - 1;
	pt::u32 E_j = r_limit_right == pc::INVALID_IDX
			      ? r_limit_right
			      : r_ref_start + r_limit_right + 1;
	pt::u32 E_i = (f_limit_right == pc::INVALID_IDX)
			      ? f_limit_right
			      : f_ref_start + f_limit_right + 1;

	// extend a
	for (; i < E_i || j > S_j; ++i, j--)
		if (do_match(f_ref_w, r_ref_w, i, j, dist_a) != 0)
			break;

	i = f_ref_start;
	j = r_ref_start;
	pt::u32 dist_b{}; // distance extended
	for (; i > S_i || j < E_j; --i, ++j)
		if (do_match(f_ref_w, r_ref_w, i, j, dist_b) != 0)
			break;

	pt::u32 x = dist_b == 0 ? f_ref_start : (f_ref_start - dist_b) + 1;
	pt::u32 y = dist_b == 0 ? r_ref_start : r_ref_start + (dist_b - 1);
	pt::u32 z = dist_a + dist_b;

	if (dist_b > 0 && dist_a > 0)
		z--;

	if (dist_a == 0 && dist_b == 0) {
		WARN("extend: no ext for ({}, {})", ref_r_id, alt_r_id);
		return std::nullopt;
	}

	return extension{f_r_id, r_r_id, x, y, z};
}

std::vector<extension> extend_links(const bd::VG &g, const chain_t &c,
				    pt::u32 ref_r_id, pt::u32 alt_r_id)
{
	std::vector<extension> extensions;
	for (pt::u32 i{}; i < c.len(); ++i) {
		const auto &lp = c.links[i];
		auto opt_ext = extend_link(g, lp, i, ref_r_id, alt_r_id);
		if (opt_ext)
			extensions.push_back(opt_ext.value());
	}

	return extensions;
}

void sne(const bd::VG &g, const std::vector<pin_cushion> &pcushions,
	 const std::set<pt::u32> &to_call_ref_ids, std::vector<pga::Exp> &exps)
{
	std::set<pt::u32> alt_ref_ids =
		find_alt_ref_ids(pcushions, to_call_ref_ids);

	// 1ref and alt
	std::set<pt::op_t<pt::u32>> call_pair;
	for (pt::u32 ref_r_id : to_call_ref_ids)
		for (pt::u32 alt_r_idx : alt_ref_ids)
			call_pair.insert({ref_r_id, alt_r_idx});

	for (auto [ref_r_id, alt_r_id] : call_pair) {
		chain_t colinear_chain =
			gen_colinear_chain(ref_r_id, alt_r_id, pcushions);

		if (colinear_chain.empty())
			continue;

		set_chain_link_limits(colinear_chain);

		std::vector<extension> e =
			extend_links(g, colinear_chain, ref_r_id, alt_r_id);

		if (e.empty())
			continue;

		std::vector<extension> merged_extensions =
			merge_extensions(e, g);

		comp_expeditions(g, merged_extensions, exps);
	}

	return;
}
} // namespace povu::overlay::sne
