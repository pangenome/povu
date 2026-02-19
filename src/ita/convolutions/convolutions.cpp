#include <cassert>
#include <thread>
#include <tuple>
#include <vector>

#include <liteseq/refs.h> // for ref_walk, ref

#include "ita/convolutions/at_matrix.hpp" // for
#include "ita/genomics/allele.hpp"	  // for hap_slice, trek
#include "ita/variation/rov.hpp"	  // for RoV
#include "meza/matrix/matrix.hpp"	  // for matrix2d
#include "povu/common/core.hpp"		  // for pt
#include "povu/graph/bidirected.hpp"	  // for VG
#include "povu/graph/types.hpp"		  // for id_or_t

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

ia::rov_boundaries
indexes_to_rov_boundaries(const bd::VG &g, const pt::op_t<pt::u32> &filter_idxs,
			  const std::vector<pt::u32> &sorted_vertices,
			  pt::u32 h_idx)
{

	const lq::ref_walk *hw = g.get_ref_vec(h_idx)->walk; // the hap walk

	auto foo = [&](pt::u32 i) -> ptg::id_or_t
	{
		pt::u32 v_id = hw->v_ids[i];
		// hw->strands[context_start];
		ptg::or_e o = pr::lq_strand_to_pv_or(hw->strands[i]);

		return {v_id, o};
	};

	auto [u, v] = filter_idxs;
	pt::u32 u_v_id = sorted_vertices[u];
	pt::u32 v_v_id = sorted_vertices[v];

	const std::vector<pt::u32> &u_positions =
		g.get_vertex_ref_idxs(g.v_id_to_idx(u_v_id), h_idx);

	const std::vector<pt::u32> &v_positions =
		g.get_vertex_ref_idxs(g.v_id_to_idx(v_v_id), h_idx);

	if (u_positions.empty() || v_positions.empty())
		std::cerr << h_idx << " u positions " << u_positions.size()
			  << " v positions " << v_positions.size() << "\n";

	pt::u32 context_start =
		std::min(u_positions.front(), v_positions.front());

	pt::u32 context_end =
		std::max(u_positions.front(), v_positions.front());

	return ia::rov_boundaries{foo(context_start), foo(context_end)};
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

	// std::vector<ia::hap_slice> alt_slices;
	for (pt::u32 h_idx : h_idxs) {
		ia::hap_slice alt_slice =
			context_to_hap_slice(g, h_idx, graph_context);
		// alt_slices.emplace_back(alt_slice);

		alt_haps.insert(alt_slice.ref_idx);

		if (alt_slice.len == ref_slice.len) // sub
			as.add_sub(std::move(alt_slice));
		else if (alt_slice.len > ref_slice.len) // ins
			as.add_ins(std::move(alt_slice));
		else // del
			as.add_del(std::move(alt_slice));
	}

	return {ref_slice, as, alt_haps};
}

std::set<pt::u32> find_ref_matches(const meza::matrix::at_matrix &result_matrix,
				   const std::set<pt::u32> &blank_rows,
				   pt::u32 I, pt::u32 ref_h_idx)
{
	auto is_match_ref = [&](pt::u32 i) -> bool
	{
		bool x = pv_cmp::contains(blank_rows, i);
		bool y = result_matrix.base().is_row_blank(i);
		// bool z = i == ref_h_idx;

		return !x && y;
	};

	std::set<pt::u32> matches_ref;
	for (pt::u32 i{}; i < I; i++)
		if (is_match_ref(i))
			matches_ref.insert(i);

	return matches_ref;
}

std::set<pt::u32>
find_variant_refs(const meza::matrix::at_matrix &result_matrix,
		  const std::set<pt::u32> &blank_rows, pt::u32 I,
		  pt::u32 ref_h_idx)
{
	auto is_variant_ref = [&](pt::u32 i) -> bool
	{
		bool x = pv_cmp::contains(blank_rows, i);
		bool y = result_matrix.base().is_row_blank(i);
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
find_context(const meza::matrix::at_matrix &ref_matrix,
	     const meza::matrix::at_matrix &result_matrix,
	     const std::set<pt::u32> &variant_refs)
{
	std::map<pt::u32, std::vector<pt::op_t<pt::u32>>> hap_contexts;
	for (pt::u32 h_idx : variant_refs)
		hap_contexts.emplace(
			h_idx, meza::matrix::find_context(ref_matrix.base(),
							  result_matrix.base(),
							  h_idx));

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

comp_res run_convolutions(const meza::matrix::at_matrix &ref_matrix,
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

	meza::matrix::matrix2d<pt::u32> res = meza::matrix::vector_xor(
		ref_matrix.base(), filter_matrix.base(), blank_rows);
	auto result_matrix = meza::matrix::at_matrix::create_from_base_matrix(
		std::move(res), filter_matrix);

	std::cerr << "R_x\n";
	result_matrix.dbg_print(std::cerr, true);

	std::set<pt::u32> matches_ref =
		find_ref_matches(result_matrix, blank_rows, I, ref_h_idx);

	std::set<pt::u32> variant_refs =
		find_variant_refs(result_matrix, blank_rows, I, ref_h_idx);

	std::map<pt::u32, std::vector<pt::op_t<pt::u32>>> hap_contexts =
		find_context(ref_matrix, result_matrix, variant_refs);

	std::map<pt::op_t<pt::u32>, std::set<pt::u32>> contexts2haps =
		contexts_to_hap_idxs(hap_contexts);

	std::map<pt::up_t<pt::u32>, std::vector<pt::u32>> at_prefix_sums =
		compare_allele_traversals(variant_refs, filter_matrix);

	// #ifdef DEBUG

	std::cerr << "---------------------------\n";
	std::cerr << "matches ref: " << pu::concat_with(matches_ref, ',')
		  << "\n";
	std::cerr << "variant refs: " << pu::concat_with(variant_refs, ',')
		  << "\n";

	std::cerr << "contexts:\n";
	for (const auto &[h_idx, contexts] : hap_contexts) {
		std::cerr << "Hap " << h_idx << "\n";
		std::cerr << "Contexts: ";
		for (auto [u, v] : contexts)
			std::cerr << "(" << u << "," << v << ")"
				  << ", ";
		std::cerr << "\n";
	}
	std::cerr << "---------------------------\n";

	// #endif

	return {contexts2haps, matches_ref, at_prefix_sums};
}

std::optional<ia::trek>
comp_expedition(const bd::VG &g, ir::RoV &rov,
		const std::set<pt::u32> &to_call_ref_ids)
{
	const pt::u32 I = g.get_hap_count();
	auto [ref_matrices, filter_matrix, result_matrix] =
		ita::at_matrix::init_depth_matrices(g, rov, to_call_ref_ids);

	const std::vector<pt::u32> &sorted_vertices = rov.get_sorted_vertices();

	auto tk = ia::trek::create_new(&rov, g.get_hap_count(), false);

	std::cerr << "R_f\n";
	filter_matrix.dbg_print(std::cerr, true);

	for (pt::u32 ref_h_idx : to_call_ref_ids) {
		const meza::matrix::depth_matrix &ref_matrix =
			ref_matrices.at(ref_h_idx);

		std::cerr << "R_r\n";
		ref_matrix.dbg_print(std::cerr, true);

		auto [context_to_haps, matches_ref, at_prefix_sums] =
			run_convolutions(ref_matrix, filter_matrix, ref_h_idx);

		for (const auto &[context, haps] : context_to_haps) {
			std::cerr << "Context: (" << context.first << ", "
				  << context.second << ")\n";

			ia::rov_boundaries ia_cxt = indexes_to_rov_boundaries(
				g, context, sorted_vertices, ref_h_idx);

			auto [ref_slice, alt_set, alt_haps] = gen_hap_slices(
				g, filter_matrix, context, ref_h_idx, haps);

			std::cerr << ref_slice.dbg_str() << "\n";

			for (auto &sub_slices : alt_set.get_subs())
				for (auto &sl : sub_slices.second)
					std::cerr << sl.dbg_str() << "\n";

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

	return tk;
}
} // namespace ita::convolutions
