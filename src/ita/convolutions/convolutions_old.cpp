#include <cassert>
#include <optional>
#include <thread>
#include <tuple>
#include <vector>

#include <liteseq/refs.h> // for ref_walk, ref

#include "ita/convolutions/at_matrix.hpp" // for matrix_pool, rov_matrix_pool, init_depth_matrices
#include "ita/genomics/allele.hpp" // for hap_slice, trek
#include "ita/variation/rov.hpp"   // for RoV
#include "meza/matrix/matrix.hpp"  // for matrix2d
#include "povu/common/compat.hpp"
#include "povu/common/constants.hpp"
#include "povu/common/core.hpp" // for pt
#include "povu/common/utils.hpp"
#include "povu/graph/bidirected.hpp" // for VG
#include "povu/graph/types.hpp"	     // for id_or_t

namespace ita::convolutions
{
namespace lq = liteseq;

struct comp_res {
	std::map<pt::op_t<pt::u32>, std::set<pt::u32>> contexts2haps;
	std::set<pt::u32> matches_ref;

	// allele traversal XORs and prefix sums
	std::map<pt::up_t<pt::u32>, std::vector<pt::u32>> at_prefix_sums;
};

pt::op_t<pt::u32>
matrix_context_to_graph_context(const meza::matrix::at_matrix &filter_matrix,
				pt::op_t<pt::u32> context)
{
	auto [u, v] = context;
	const std::vector<pt::u32> &col_names = filter_matrix.get_col_names();
	pt::u32 u_v_id = col_names[u];
	pt::u32 v_v_id = col_names[v];

	return {u_v_id, v_v_id};
}

std::map<pt::op_t<pt::u32>, std::set<pt::u32>> contexts_to_hap_idxs(
	const std::map<pt::u32, std::vector<pt::op_t<pt::u32>>> &hap_contexts)
{
	std::map<pt::op_t<pt::u32>, std::set<pt::u32>> res;
	for (const auto &[h_idx, contexts] : hap_contexts)
		for (pt::op_t<pt::u32> context : contexts)
			res[context].insert(h_idx);

	return res;
}

ia::rov_boundaries indexes_to_rov_boundaries_no_tangle(
	const bd::VG &g, const pt::op_t<pt::u32> &filter_idxs,
	const std::vector<pt::u32> &sorted_vertices, pt::u32 h_idx)
{
	const lq::ref_walk *hw = g.get_ref_vec(h_idx)->walk; // the hap walk

	auto step_idx2id_or = [&](pt::u32 i) -> ptg::id_or_t
	{
		pt::u32 v_id = hw->v_ids[i];
		ptg::or_e o = pr::lq_strand_to_pv_or(hw->strands[i]);
		return {v_id, o};
	};

	auto [u, v] = filter_idxs;
	pt::u32 u_v_idx = g.v_id_to_idx(sorted_vertices[u]);
	pt::u32 v_v_idx = g.v_id_to_idx(sorted_vertices[v]);

	pt::u32 u_step_idx = g.get_vertex_ref_idxs(u_v_idx, h_idx).front();
	pt::u32 v_step_idx = g.get_vertex_ref_idxs(v_v_idx, h_idx).back();

	ptg::id_or_t u_id_or = step_idx2id_or(u_step_idx);
	ptg::id_or_t v_id_or = step_idx2id_or(v_step_idx);

	return ia::rov_boundaries{u_id_or, v_id_or};
}

ia::rov_boundaries indexes_to_rov_boundaries_tangled(
	const bd::VG &g, const pt::op_t<pt::u32> &filter_idxs,
	const std::vector<pt::u32> &sorted_vertices,
	const ita::traversals::traversals::allele_traversal &at, pt::u32 h_idx)
{
	const lq::ref_walk *hw = g.get_ref_vec(h_idx)->walk; // the hap walk

	auto step_idx2id_or = [&](pt::u32 i) -> ptg::id_or_t
	{
		pt::u32 v_id = hw->v_ids[i];
		ptg::or_e o = pr::lq_strand_to_pv_or(hw->strands[i]);
		return {v_id, o};
	};

	auto [u, v] = filter_idxs;
	pt::u32 u_v_id = sorted_vertices[u];
	pt::u32 v_v_id = sorted_vertices[v];

	pt::u32 u_step_idx = pc::MAX_IDX;
	pt::u32 v_step_idx = 0;

	if (at.empty())
		throw std::runtime_error("allele traversal is empty");

	for (pt::u32 step_idx : at) {
		pt::u32 v_id = hw->v_ids[step_idx];

		if (step_idx < u_step_idx && v_id == u_v_id)
			u_step_idx = step_idx;

		if (v_step_idx < step_idx && v_id == v_v_id)
			v_step_idx = step_idx;
	}

	ptg::id_or_t u_id_or = step_idx2id_or(u_step_idx);
	ptg::id_or_t v_id_or = step_idx2id_or(v_step_idx);

	return ia::rov_boundaries{u_id_or, v_id_or};
}

ia::hap_slice context_to_hap_slice(const bd::VG &g, pt::u32 h_idx,
				   pt::op_t<pt::u32> graph_context)
{
	auto [u_v_id, v_v_id] = graph_context;
	const std::vector<pt::u32> &u_positions =
		g.get_vertex_ref_idxs(g.v_id_to_idx(u_v_id), h_idx);

	const std::vector<pt::u32> &v_positions =
		g.get_vertex_ref_idxs(g.v_id_to_idx(v_v_id), h_idx);

#ifdef DEBUG
	assert(u_positions.size() == v_positions.size());
	assert(!u_positions.empty());
	assert(!v_positions.empty());
#endif

	pt::u32 slice_start =
		std::min(u_positions.front(), v_positions.front());
	pt::u32 b = std::max(u_positions.front(), v_positions.front());
	pt::u32 len = (b - slice_start) + 1;

	return {g.get_ref_vec(h_idx)->walk, h_idx, slice_start, len};
}

std::tuple<ia::hap_slice, ia::alt_set, std::set<pt::u32>>
gen_hap_slices(const bd::VG &g, const meza::matrix::at_matrix &filter_matrix,
	       pt::op_t<pt::u32> context, pt::u32 ref_h_idx,
	       const std::set<pt::u32> &h_idxs)
{
	pt::op_t<pt::u32> graph_context =
		matrix_context_to_graph_context(filter_matrix, context);

	ia::hap_slice ref_slice =
		context_to_hap_slice(g, ref_h_idx, graph_context);

	ia::alt_set as; // alt set

	std::set<pt::u32> alt_haps;

	for (pt::u32 h_idx : h_idxs) {
		ia::hap_slice alt_slice =
			context_to_hap_slice(g, h_idx, graph_context);

		alt_haps.insert(alt_slice.ref_idx);

		if (alt_slice.len == ref_slice.len) // ....... sub
			as.add_sub(std::move(alt_slice));
		else if (alt_slice.len > ref_slice.len) // ... ins
			as.add_ins(std::move(alt_slice));
		else // ...................................... del
			as.add_del(std::move(alt_slice));
	}

	return {ref_slice, as, alt_haps};
}

/**
 * @brief includes the ref hap idx
 */
std::set<pt::u32>
find_ref_matches(const meza::matrix::dense_matrix2d<pt::u32> &result_matrix,
		 const std::set<pt::u32> &blank_rows, pt::u32 I)
{
	auto is_match_ref = [&](pt::u32 i) -> bool
	{
		bool x = pv_cmp::contains(blank_rows, i);
		bool y = result_matrix.is_row_blank(i);

		return !x && y;
	};

	std::set<pt::u32> matches_ref;
	for (pt::u32 i{}; i < I; i++)
		if (is_match_ref(i))
			matches_ref.insert(i);

	return matches_ref;
}

std::set<pt::u32>
find_variant_refs(const meza::matrix::dense_matrix2d<pt::u32> &result_matrix,
		  const std::set<pt::u32> &blank_rows, pt::u32 I,
		  pt::u32 ref_h_idx)
{
	auto is_variant_ref = [&](pt::u32 i) -> bool
	{
		bool x = pv_cmp::contains(blank_rows, i);
		bool y = result_matrix.is_row_blank(i);
		bool z = i == ref_h_idx;

		return !x && !y && !z;
	};

	std::set<pt::u32> variant_refs;
	for (pt::u32 i{}; i < I; i++)
		if (is_variant_ref(i))
			variant_refs.insert(i);

	return variant_refs;
}

std::map<pt::u32, std::vector<pt::op_t<pt::u32>>>
find_context2(const meza::matrix::ref_matrix &ref_matrix,
	      const meza::matrix::dense_matrix2d<pt::u32> &result_matrix,
	      const std::set<pt::u32> &variant_refs)
{
	std::map<pt::u32, std::vector<pt::op_t<pt::u32>>> hap_contexts;
	for (pt::u32 h_idx : variant_refs) {
		std::vector<pt::op_t<pt::u32>> contexts =
			meza::matrix::find_context(ref_matrix.base(),
						   result_matrix, h_idx);
		hap_contexts.emplace(h_idx, std::move(contexts));
	}

	return hap_contexts;
}

std::map<pt::up_t<pt::u32>, std::vector<pt::u32>>
compare_allele_traversals(const std::set<pt::u32> &variant_refs,
			  const meza::matrix::at_matrix &filter_matrix)
{
	std::map<pt::up_t<pt::u32>, std::vector<pt::u32>> m;
	for (pt::u32 h1_idx : variant_refs) {
		for (pt::u32 h2_idx : variant_refs) {
			pt::up_t<pt::u32> k{h1_idx, h2_idx}; // pair key
			if (h1_idx == h2_idx || pv_cmp::contains(m, k))
				continue;

			// prefix sum of xor row, for allele traversal
			std::vector<pt::u32> ps =
				filter_matrix.base().xor_rows_and_prefix_sum(
					h1_idx, h2_idx);

			m[k] = ps;
		}
	}

	return m;
}

comp_res run_convolutions(const meza::matrix::ref_matrix &ref_matrix,
			  const meza::matrix::at_matrix &filter_matrix,
			  pt::u32 ref_h_idx)
{
	auto is_row_blank = [&](pt::u32 i)
	{
		return filter_matrix.base().is_row_blank(i);
	};

	std::set<pt::u32> blank_rows;
	pt::u32 I = ref_matrix.base().rows();
	for (pt::u32 i{}; i < I; i++)
		if (is_row_blank(i))
			blank_rows.insert(i);

	// std::cerr << "ref\n";
	// ref_matrix.base().dbg_print(std::cerr);

	// std::cerr << "filter\n";
	// filter_matrix.base().dbg_print(std::cerr);

	meza::matrix::dense_matrix2d<pt::u32> res = meza::matrix::vector_xor(
		ref_matrix.base(), filter_matrix.base(), blank_rows);

	// std::cerr << "res\n";
	// res.dbg_print(std::cerr);

	std::set<pt::u32> matches_ref = find_ref_matches(res, blank_rows, I);

	std::set<pt::u32> variant_refs =
		find_variant_refs(res, blank_rows, I, ref_h_idx);

	std::map<pt::u32, std::vector<pt::op_t<pt::u32>>> hap_contexts =
		find_context2(ref_matrix, res, variant_refs);

	std::map<pt::op_t<pt::u32>, std::set<pt::u32>> contexts2haps =
		contexts_to_hap_idxs(hap_contexts);

	std::map<pt::up_t<pt::u32>, std::vector<pt::u32>> at_prefix_sums =
		compare_allele_traversals(variant_refs, filter_matrix);

	// std::cerr << "match ref: " << pu::concat_with(matches_ref, ',') <<
	// "\n"; std::cerr << "var ref: " << pu::concat_with(variant_refs, ',')
	// << "\n";

	return {contexts2haps, matches_ref, at_prefix_sums};
}

void conv_pool(
	const bd::VG &g, const std::vector<pt::u32> &sorted_vertices,
	const std::set<pt::u32> &to_call_ref_ids,
	ita::at_matrix::matrix_pool &pool,
	const std::vector<ita::traversals::traversals::itinerary> &hap_itns,
	ia::trek &tk)
{
	const auto &[ref_matrices, filter_matrix, _, hap_slices, ref_loop_no,
		     is_tangled, I, __] = pool;

	for (pt::u32 ref_h_idx : to_call_ref_ids) {

		if (!pv_cmp::contains(ref_matrices, ref_h_idx))
			continue;

		const meza::matrix::ref_matrix &ref_matrix =
			ref_matrices.at(ref_h_idx);

		auto [context_to_haps, matches_ref, at_prefix_sums] =
			run_convolutions(ref_matrix, filter_matrix, ref_h_idx);

		for (const auto &[context, haps] : context_to_haps) {

			std::optional<ia::rov_boundaries> opt_ia_cxt;
			if (is_tangled) {
				const ita::traversals::traversals::itinerary
					&ref_hap_itn = hap_itns.at(ref_h_idx);
				const ita::traversals::traversals::
					allele_traversal &at =
						ref_hap_itn.at(ref_loop_no);
				opt_ia_cxt = indexes_to_rov_boundaries_tangled(
					g, context, sorted_vertices, at,
					ref_h_idx);
			}
			else {
				opt_ia_cxt =
					indexes_to_rov_boundaries_no_tangle(
						g, context, sorted_vertices,
						ref_h_idx);
			}

			ia::rov_boundaries ia_cxt = opt_ia_cxt.value();

			auto [ref_slice, alt_set, alt_haps] = gen_hap_slices(
				g, filter_matrix, context, ref_h_idx, haps);

			ia::cxt_to_min_rov_map d;
			d.emplace(ia_cxt,
				  ia::minimal_rov(ia_cxt, ref_slice,
						  std::move(matches_ref),
						  std::move(alt_haps),
						  std::move(alt_set)));

			tk.set_min_rov(ref_h_idx, std::move(d));

			for (auto h_idx : matches_ref)
				tk.add_match_ref(ref_h_idx, h_idx);

			for (pt::u32 h_idx{}; h_idx < I; h_idx++)
				if (filter_matrix.base().is_row_blank(h_idx) &&
				    h_idx != ref_h_idx)
					tk.add_no_cov(ref_h_idx, h_idx);
		}
	}
}

std::optional<ia::trek>
comp_expedition(const bd::VG &g, ir::RoV &rov,
		const std::set<pt::u32> &to_call_ref_ids)
{
	using ita::at_matrix::matrix_pool;
	using ita::at_matrix::rov_matrix_pool;

	rov_matrix_pool rov_mp =
		ita::at_matrix::init_depth_matrices(g, rov, to_call_ref_ids);

	if (rov_mp.empty())
		return std::nullopt;

	auto &[_, pools, is_tangled, hap_itns] = rov_mp;

	const std::vector<pt::u32> &sorted_vertices = rov.get_sorted_vertices();

	auto tk = ia::trek::create_new(&rov, g.get_hap_count(), is_tangled);

	if (is_tangled) {
		for (matrix_pool &pool : pools)
			conv_pool(g, sorted_vertices, to_call_ref_ids, pool,
				  hap_itns, tk);
	}
	else { //  always only one pool if not tangled
		matrix_pool &pool = pools.front();
		conv_pool(g, sorted_vertices, to_call_ref_ids, pool, hap_itns,
			  tk);
	}

	if (!tk.has_data())
		return std::nullopt;

	return tk;
}
} // namespace ita::convolutions
