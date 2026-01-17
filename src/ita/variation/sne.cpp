#include "ita/variation/sne.hpp"

// std includes
#include <algorithm> // for min, max
#include <optional>

// deps
#include <liteseq/refs.h> // for ref_walk, ref

// ita includes
#include "ita/graph/slice_tree.hpp" // for ist

// povu includes
#include "povu/common/constants.hpp"
#include "povu/common/core.hpp"
#include "povu/graph/pvst.hpp"
#include "povu/graph/types.hpp" // pgt

namespace ita::sne
{

namespace lq = liteseq;
namespace pgt = povu::types::graph;

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
 * Adds 1 before division for ceiling behavior
 */
pt::u32 half_diff_ceil(pt::u32 a, pt::u32 b)
{
	pt::u32 diff = std::max(a, b) - std::min(a, b);
	return (diff + 1) / 2;
}

std::optional<chain_t> gen_colinear_chain(pt::u32 ref_h_idx, pt::u32 alt_h_idx,
					  const pin_cushion &pc)
{
	chain_t chain;

	pt::up_t<pt::u32> k{ref_h_idx, alt_h_idx};
	std::vector<std::pair<pin, pin>> y = pc.get_pin_pairs(k);

	if (y.empty())
		return std::nullopt;

	for (auto [a_pin, b_pin] : y) {
		auto [ref_pin, alt_pin] =
			(a_pin.r_idx == ref_h_idx)
				? std::make_pair(a_pin, b_pin)
				: std::make_pair(b_pin, a_pin);

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
				     pt::u32 ref_r_id, pt::u32 alt_r_id)
{
	auto [ref_link, alt_link] = lp;

	// auto [fwd_link, rev_link] =
	//	(ref_link.dir == hap_traversal::left)
	//		? std::make_pair(ref_link, alt_link)
	//		: std::make_pair(alt_link, ref_link);

	// auto [f_r_id, r_r_id] = (ref_link.dir == hap_traversal::left)
	//				? std::make_pair(ref_r_id, alt_r_id)
	//				: std::make_pair(alt_r_id, ref_r_id);

	const lq::ref_walk *ref_h_w = g.get_ref_vec(ref_r_id)->walk;
	const lq::ref_walk *alt_h_w = g.get_ref_vec(alt_r_id)->walk;

	// INFO("hee");
	// std::cerr << "(ref: " << ref_h_w->v_ids[ref_link.idx] << ", ";
	// std::cerr << "alt: " << alt_h_w->v_ids[alt_link.idx] << ")\n";

	// const lq::ref_walk *f_ref_w = g.get_ref_vec(f_r_id)->walk;
	// const lq::ref_walk *r_ref_w = g.get_ref_vec(r_r_id)->walk;

	// auto [_, f_limit_left, f_limit_right, _] = fwd_link;
	// auto [_, r_limit_left, r_limit_right, __] = rev_link;

	pt::u32 f_limit_left = ref_link.limit_left;
	pt::u32 f_limit_right = ref_link.limit_right;

	pt::u32 r_limit_left = alt_link.limit_left;
	pt::u32 r_limit_right = alt_link.limit_right;

	pt::u32 f_ref_start = ref_link.idx;
	pt::u32 r_ref_start = alt_link.idx;

	pt::u32 i{ref_link.idx};
	pt::u32 j{alt_link.idx};
	pt::u32 dist_a{}; // distance extended

	// i is ref; j is alt
	pt::u32 S_i = (f_limit_left == pc::INVALID_IDX)
			      ? f_ref_start
			      : f_ref_start - f_limit_left - 1;
	pt::u32 S_j = (f_limit_right == pc::INVALID_IDX)
			      ? r_ref_start
			      : r_ref_start - r_limit_left - 1;
	pt::u32 E_j = r_limit_right == pc::INVALID_IDX
			      ? r_limit_right
			      : r_ref_start + r_limit_right + 1;
	pt::u32 E_i = (f_limit_right == pc::INVALID_IDX)
			      ? f_limit_right
			      : f_ref_start + f_limit_right + 1;

	// std::cerr << "S_i" << S_i << " r " << f_ref_start << " E_i " << E_i
	//	  << "\n"
	//	  << "S_j" << S_j << " a " << r_ref_start << " E_j " << E_j
	//	  << "\n";

	// std::cerr << "S_i" << S_i << "S_j" << S_j << "E_i" << E_i << "E_j"
	//	  << E_j << "\n";

	// extend a
	for (; i < E_i || j > S_j; ++i, j--)
		if (do_match(ref_h_w, alt_h_w, i, j, dist_a) != 0)
			break;

	i = f_ref_start;
	j = r_ref_start;
	pt::u32 dist_b{}; // distance extended
	for (; i > S_i || j < E_j; --i, ++j)
		if (do_match(ref_h_w, alt_h_w, i, j, dist_b) != 0)
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

	// std::cerr << "Dist a" << dist_a << "dist b" << dist_b << "\n";

	return extension{ref_r_id, alt_r_id, x, y, z};
}

void extend_links(const bd::VG &g, const chain_t &c, pt::u32 ref_r_id,
		  pt::u32 alt_r_id, ist::st &it_)
{
	// bool dbg = (alt_r_id == 83) ? true : false;
	// dbg = false;

	std::vector<extension> extensions;
	for (pt::u32 i{}; i < c.len(); ++i) {
		const auto &lp = c.links[i];
		auto opt_ext = extend_link(g, lp, ref_r_id, alt_r_id);
		if (opt_ext) {
			extensions.push_back(opt_ext.value());
			auto [_, r_r_idx, ref_h_start, r_ref_start, len] =
				opt_ext.value();

			pt::u32 alt_h_start = r_ref_start - len + 1;

			// if (true) {
			//	std::cerr << "\n";
			//	std::cerr << "ref start  " << ref_h_start
			//		  << ", alt hap " << r_r_idx
			//		  << ", alt start " << alt_h_start
			//		  << ", len " << len << "\n";
			//	opt_ext.value().dbg_print(g);
			// }

			// for (const auto &[id, v] :
			// it_.get_vertices()) {
			//	if (!v.has_alts()) {
			//		ERR("Vertex {} has no alts",
			// id);		std::exit(EXIT_FAILURE);
			//	}
			// }

			// if (dbg) {
			//	std::cerr << "ADDING TO ITREE:
			// ref_h_start "
			//		  << ref_h_start << ",
			// alt_h_start "
			//		  << alt_h_start << ", len " <<
			// len
			//		  << "\n";

			//	std::cerr << "ITREE BEFORE ADDING:\n";
			//	it_.dbg_print(std::cerr, g);
			// }

			it_.add_vertex(ref_h_start, r_r_idx, alt_h_start, len);

			// if (pv_cmp::contains(it_.get_vertices(), 63)) {
			//	const ist::vertex &v =
			//		it_.get_vertices().at(63);
			//	INFO("Vertex 63 exists in itree");
			//	if (!v.has_alts()) {
			//		ERR("Vertex 63 has no alts");
			//		std::exit(EXIT_FAILURE);
			//	}
			//	else {
			//		for (auto &[len, alts] :
			//		     v.get_same_len_alts()) {
			//			INFO("Vertex 63 alt len {} "
			//			     "count {}",
			//			     len, alts.size());
			//		}
			//	}
			// }

			// if (dbg) {
			//	std::cerr << "ITREE AFTER ADDING:\n";
			//	it_.dbg_print(std::cerr, g);
			// }
		}
	}

	// std::exit(0);

	// INFO("Extended {} links", extensions.size());
	// for (const auto &ext : extensions) {
	//	INFO("EXT");
	//	ext.dbg_print(g);
	// }

	// return;
}

std::vector<pt::slice> find_hap_slices(const bd::VG &g, pt::u32 h_idx,
				       pt::u32 u_v_id, pt::u32 v_v_id)
{
	std::vector<pt::u32> u_positions =
		g.get_vertex_ref_idxs(g.v_id_to_idx(u_v_id), h_idx);

	std::vector<pt::u32> v_positions =
		g.get_vertex_ref_idxs(g.v_id_to_idx(v_v_id), h_idx);

	if (u_positions.empty() || v_positions.empty())
		return {};

	std::sort(u_positions.begin(), u_positions.end());
	std::sort(v_positions.begin(), v_positions.end());

	std::vector<pt::slice> slices;
	pt::u32 N = std::min(u_positions.size(), v_positions.size());

	for (pt::u32 i{}; i < N; ++i) {
		pt::u32 x = std::min(u_positions[i], v_positions[i]);
		pt::u32 y = std::max(u_positions[i], v_positions[i]);

		if (x == y)
			continue;

		slices.emplace_back(x, y - x + 1);
	}

	return slices;
}

bool is_slice_in_hap(const bd::VG &g, pt::u32 ref_h_idx, pt::u32 ref_h_start,
		     pt::u32 len, pt::u32 h_idx)
{
	const lq::ref_walk *h_w1 = g.get_ref_vec(ref_h_idx)->walk;
	const lq::ref_walk *h_w2 = g.get_ref_vec(h_idx)->walk;

	pt::u32 s = ref_h_start;
	pt::u32 N = ref_h_start + len;

	pt::u32 u_v_id = h_w1->v_ids[s];
	pt::u32 v_v_id = h_w1->v_ids[s];

	std::vector<pt::slice> h2_slices =
		find_hap_slices(g, h_idx, u_v_id, v_v_id);

	if (h2_slices.empty())
		return false;

	auto get_at = [](const lq::ref_walk *h_w, pt::u32 i) -> bd::id_or_t
	{
		pt::idx_t ref_v_id = h_w->v_ids[i];
		pgt::or_e ref_o = h_w->strands[i] == lq::strand::STRAND_FWD
					  ? pgt::or_e::forward
					  : pgt::or_e::reverse;
		bd::id_or_t step{ref_v_id, ref_o};
		return step;
	};

	for (auto sl : h2_slices) {
		if (sl.len() != len)
			continue;

		pt::u32 j = sl.start();

		for (pt::u32 i{s}; i < N; i++, j++)
			if (get_at(h_w1, i) != get_at(h_w2, j))
				break;

		if (j == sl.start() + sl.len() - 1)
			return true; // we got to the end, matched
	}

	return false;
}

// void find_ref_haps(const bd::VG &g, ist::st &i_tree, pt::u32 I)
// {
//	// pt::u32 N = i_tree.size();
//	pt::u32 ref_h_idx = i_tree.get_ref_hap_idx();

//	for (const auto &[_, v] : i_tree.get_vertices()) {
//		for (auto &[len, alts] : v.get_same_len_alts()) {
//			// INFO("LEN {}", len);

//			v.add_len_haps(len, ref_h_idx);

//			std::set<pt::u32> alt_haps;
//			for (const auto &a : alts)
//				alt_haps.insert(a.h_idx);

//			for (pt::u32 h_idx{}; h_idx < I; ++h_idx) {
//				if (pv_cmp::contains(alt_haps, h_idx))
//					continue;

//				if (h_idx == ref_h_idx)
//					continue;

//				if (is_slice_in_hap(g, ref_h_idx,
// v.ref_h_start,						    len,
// h_idx)) {					v.add_len_haps(h_idx,
// len);
//				}
//			}
//		}
//	}

//	// for (pt::u32 it_v_idx{}; it_v_idx < N; it_v_idx++) {
//	//	ist::vertex &v = i_tree.get_vertex_mut(it_v_idx);
//	//	// std::set<pt::u32> alt_haps = v.get_alt_hap_indices();

//	// }
// }

std::vector<ist::st> sne(const bd::VG &g, const pin_cushion &pcushions,
			 const std::set<pt::u32> &to_call_ref_ids)
{
	const pt::u32 I = g.get_hap_count(); // rows
	std::vector<ist::st> i_trees;

	for (pt::u32 ref_h_idx : to_call_ref_ids) {
		ist::st i_tree(ref_h_idx);
		for (pt::u32 alt_h_idx{}; alt_h_idx < I; ++alt_h_idx) {

			if (ref_h_idx == alt_h_idx)
				continue;

			// const povu::refs::Ref &r =
			// g.get_ref_by_id(alt_h_idx); std::cerr <<
			// "Tag: " << r.tag() << "\n";

			std::optional<chain_t> opt_co_chain =
				gen_colinear_chain(ref_h_idx, alt_h_idx,
						   pcushions);

			if (!opt_co_chain)
				continue;

			chain_t colinear_chain = *opt_co_chain;

			if (colinear_chain.empty())
				continue;

			set_chain_link_limits(colinear_chain);

			extend_links(g, colinear_chain, ref_h_idx, alt_h_idx,
				     i_tree);
		}

		if (!i_tree.is_empty())
			i_trees.emplace_back(std::move(i_tree));
	}

	return i_trees;
}
} // namespace ita::sne
