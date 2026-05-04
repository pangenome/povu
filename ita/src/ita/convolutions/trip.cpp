#include "ita/convolutions/trip.hpp"

#include <algorithm> // for std::min, std::max
#include <optional>  // for std::optional
#include <set>	     // for std::set
#include <vector>    // for std::vector

#include <liteseq/refs.h> // for ref_walk, ref
#include <log.h>
#include <meza/pool/hap_comp.hpp>	  // for haps_comp_set
#include <meza/pool/split_pool_types.hpp> // for ov_mat_t
#include <quilt/types.hpp>		  // for qt

#include "ita/traversals/at_matrix.hpp" // for matrix_pool, rov_matrix_set

// #include "quilt/types.hpp"

namespace ita::trip
{
namespace lq = liteseq;
using meza::pool::ov_mat_t;
using mat_context_t = qt::op_t<qt::u32>;
using graph_context_t = qt::op_t<qt::u32>;
using meza::pool::hap_comp::haps_comp_set;

// map of context to the set of hap indexes with that context as a variant
using cxt_idx_t = std::map<mat_context_t, std::set<qt::u32>>;

graph_context_t
matrix_context_to_graph_context(const std::vector<qt::u32> &sorted_vertices,
				const qt::op_t<qt::u32> &matrix_context)
{
	auto [u, v] = matrix_context;
	qt::u32 u_v_id = sorted_vertices[u];
	qt::u32 v_v_id = sorted_vertices[v];

	return {u_v_id, v_v_id};
}

// use when no tangling
qt::op_t<qt::u32>
find_hap_slice_no_tangle(const bd::VG &g, qt::u32 h_idx,
			 const qt::op_t<qt::u32> &graph_context)
{
	auto [u_v_id, v_v_id] = graph_context;

	qt::u32 u_step_idx{pc::INVALID_IDX};
	qt::u32 v_step_idx{pc::INVALID_IDX};

	const std::vector<qt::u32> &u_positions =
		g.get_vertex_ref_idxs(g.v_id_to_idx(u_v_id), h_idx);

	const std::vector<qt::u32> &v_positions =
		g.get_vertex_ref_idxs(g.v_id_to_idx(v_v_id), h_idx);

	if (u_positions.empty() || v_positions.empty())
		throw std::runtime_error("One of the vertices in the context "
					 "is not present in the hap walk");

	if (u_positions.size() > 1 || v_positions.size() > 1)
		log_warn(
			"Multiple positions found for a vertex in the context. "
			"This should not happen when no tangling is present.");

	u_step_idx = std::min(u_positions.front(), v_positions.front());
	v_step_idx = std::max(u_positions.front(), v_positions.front());

	return {u_step_idx, v_step_idx};
}

qt::op_t<qt::u32> find_hap_slice_tangled(const bd::VG &g, qt::u32 h_idx,
					 const std::vector<qt::u32> &at,
					 const qt::op_t<qt::u32> &graph_context)
{
	auto [u_v_id, v_v_id] = graph_context;
	const lq::ref_walk *hw = g.get_ref_vec(h_idx)->walk; // the hap walk

	qt::u32 u_step_idx{pc::INVALID_IDX};
	qt::u32 v_step_idx{pc::INVALID_IDX};

	for (qt::u32 step_idx : at) {
		if (hw->v_ids[step_idx] == u_v_id)
			u_step_idx = step_idx;

		if (hw->v_ids[step_idx] == v_v_id)
			v_step_idx = step_idx;
	}

	return {u_step_idx, v_step_idx};
}

// TODO: make this actual slice?
qt::op_t<qt::u32> find_hap_slice(
	const bd::VG &g, qt::u32 h_idx, qt::u32 loop_no,
	const std::vector<ita::traversals::traversals::itinerary> &hap_itns,
	const graph_context_t &graph_context)
{
	if (!hap_itns.empty()) {
		const ita::traversals::traversals::itinerary &ref_hap_itn =
			hap_itns.at(h_idx);
		const ita::traversals::traversals::allele_traversal &at =
			ref_hap_itn.at(loop_no);
		return find_hap_slice_tangled(g, h_idx, at, graph_context);
	}

	return find_hap_slice_no_tangle(g, h_idx, graph_context);
}

/**
 * Finds contexts in the given row of the matrix that could potentially
 * contain a variant. A context is defined as a pair of column indices
 * (u, v) such that:
 * 1. The reference matrix has a non-zero value at (i, u) and (i, v).
 * 2. The XOR matrix has a zero value at (i, u) and (i, v).
 * 3. There are at least one difference between columns u and v in the XOR
 * matrix.
 * 4. Columns u and v are not adjacent.
 *
 *
 * @param ref_matrix The reference matrix.
 * @param xor_matrix The XOR matrix.
 * @param i The row index to analyze.
 * @return A vector of contexts represented as pairs of column indices.
 */
std::vector<mat_context_t> find_row_context(const ov_mat_t &ref_matrix,
					    const ov_mat_t &xor_mat, qt::u32 i)
{
	qt::u32 J = ref_matrix.cols();
	if (J == 1)
		return {};

	std::vector<qt::u32> sum_row(J, 0);
	sum_row[0] = xor_mat.base().at(i, 0);
	for (qt::u32 j{1}; j < J; j++) {
		qt::u8 val = xor_mat.base().at(i, j);
		sum_row[j] = sum_row[j - 1] + val;
	}

	std::vector<mat_context_t> bounds;
	std::queue<qt::u32> q; // a buffer to store similar cols

	qt::u32 j{}, u{}, v{};

	for (; j < J; j++) {
		qt::u8 ref_val = ref_matrix.base().at(i, j);
		qt::u8 xor_val = xor_mat.base().at(i, j);

		if (ref_val != 0 && xor_val == 0)
			q.push(j);

		if (q.size() > 1) {
			u = q.front();
			q.pop();
			v = q.front();

			// if (u,v) are not adjacent & contain a
			// difference
			if ((v - u) > 1 && sum_row[v] - sum_row[u] > 0)
				bounds.emplace_back(u, v);
		}
	}

	return bounds;
}

// TODO: combine with find_row_context
cxt_idx_t find_variant_mat_contexts(const ita::at_matrix::mat3 &mat_set,
				    const std::set<qt::u32> &variant_refs)
{
	const ov_mat_t &ref_mat = mat_set.ref;
	const ov_mat_t &xor_mat = mat_set.xor_result;

	cxt_idx_t index;
	// variant_context_map_t hap_contexts;
	for (qt::u32 h_idx : variant_refs) {
		std::vector<mat_context_t> contexts =
			find_row_context(ref_mat, xor_mat, h_idx);
		for (mat_context_t &c : contexts)
			index[c].insert(h_idx);
		// index.add_mat_context(h_idx, std::move(c));
	}

	return index;
}

ia::rov_boundaries indexes_to_rov_boundaries(
	const bd::VG &g, const qt::op_t<qt::u32> &matrix_context,
	const std::vector<qt::u32> &sorted_vertices,
	const std::vector<ita::traversals::traversals::itinerary> &hap_itns,
	qt::u32 h_idx, const ita::at_matrix::hap2loop &h2l)
{
	const lq::ref_walk *hw = g.get_ref_vec(h_idx)->walk; // the hap walk

	graph_context_t graph_context = matrix_context_to_graph_context(
		sorted_vertices, matrix_context);

	auto [u_step_idx, v_step_idx] = find_hap_slice(
		g, h_idx, h2l.get_loop_no(h_idx), hap_itns, graph_context);

	if (u_step_idx == pc::INVALID_IDX || v_step_idx == pc::INVALID_IDX)
		throw std::runtime_error(
			"Could not find step idx for u or v in tangle");

	auto [u_v_id, v_v_id] = graph_context;
	ptg::or_e u_o = pr::lq_strand_to_pv_or(hw->strands[u_step_idx]);
	ptg::or_e v_o = pr::lq_strand_to_pv_or(hw->strands[u_step_idx]);

	ptg::id_or_t u_id_or = {u_v_id, u_o};
	ptg::id_or_t v_id_or = {v_v_id, v_o};

	return ia::rov_boundaries{u_id_or, v_id_or};
}

ia::hap_slice context_to_hap_slice(
	const bd::VG &g, qt::u32 h_idx, qt::u32 loop_no,
	const qt::op_t<qt::u32> &graph_context,
	const std::vector<ita::traversals::traversals::itinerary> &hap_itns)
{
	auto [slice_start, slice_end] =
		find_hap_slice(g, h_idx, loop_no, hap_itns, graph_context);

	qt::u32 len = (slice_end - slice_start) + 1;

	return {g.get_ref_vec(h_idx)->walk, h_idx, slice_start, len};
}

std::tuple<ia::hap_slice, ia::alt_set, std::set<qt::u32>> gen_hap_slices(
	const bd::VG &g, const std::vector<qt::u32> &sorted_vertices,
	const qt::op_t<qt::u32> &matrix_context, qt::u32 ref_h_idx,
	const ita::at_matrix::hap2loop &h2l,
	const std::vector<ita::traversals::traversals::itinerary> &hap_itns,
	const std::set<qt::u32> &h_idxs)
{
	qt::op_t<qt::u32> graph_context = matrix_context_to_graph_context(
		sorted_vertices, matrix_context);

	ia::hap_slice ref_slice =
		context_to_hap_slice(g, ref_h_idx, h2l.get_loop_no(ref_h_idx),
				     graph_context, hap_itns);

	ia::alt_set as; // alt set

	std::set<qt::u32> alt_haps;

	for (qt::u32 h_idx : h_idxs) {
		qt::u32 loop_no = h2l.get_loop_no(h_idx);
		ia::hap_slice alt_slice = context_to_hap_slice(
			g, h_idx, loop_no, graph_context, hap_itns);

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

void gen_trip(
	const bd::VG &g, const cxt_idx_t &context_idx, qt::u32 ref_h_idx,
	const ita::at_matrix::hap2loop &h2l,
	const std::vector<ita::traversals::traversals::itinerary> &hap_itns,
	const std::vector<qt::u32> &sorted_vertices,
	const std::set<qt::u32> &matches_ref, const ov_mat_t &filter_mat,
	ia::trek &tk)
{
	const qt::u32 I = g.get_hap_count();

	for (const auto &[matrix_cxt, haps] : context_idx) {
		ia::rov_boundaries ia_cxt = indexes_to_rov_boundaries(
			g, matrix_cxt, sorted_vertices, hap_itns, ref_h_idx,
			h2l);
		// ia::rov_boundaries ia_cxt =
		//	(is_tangled) ? indexes_to_rov_boundaries_tangled(
		//			       g, context, sorted_vertices,
		//			       hap_itns, ref_h_idx, h2l)
		//		     : indexes_to_rov_boundaries_no_tangle(
		//			       g, context, sorted_vertices,
		//			       ref_h_idx);

		auto [ref_slice, alt_set, alt_haps] =
			gen_hap_slices(g, sorted_vertices, matrix_cxt,
				       ref_h_idx, h2l, hap_itns, haps);

		ia::cxt_to_min_rov_map d;
		d.emplace(ia_cxt,
			  ia::minimal_rov(ia_cxt, ref_slice, matches_ref,
					  alt_haps, std::move(alt_set)));

		tk.set_min_rov(ref_h_idx, std::move(d));

		for (auto h_idx : matches_ref)
			tk.add_match_ref(ref_h_idx, h_idx);

		for (qt::u32 h_idx{}; h_idx < I; h_idx++)
			if (filter_mat.base().is_row_blank(h_idx) &&
			    h_idx != ref_h_idx)
				tk.add_no_cov(ref_h_idx, h_idx);
	}
}

std::optional<ia::trek>
gen_trip(const bd::VG &g, const ir::RoV *rov, bool is_tangled,
	 qt::u32 ref_h_idx, const ita::at_matrix::hap2loop &h2l,
	 const std::vector<ita::traversals::traversals::itinerary> &hap_itns,
	 const ita::at_matrix::mat3 &mat_set,
	 const std::vector<qt::u32> &sorted_vertices,
	 const haps_comp_set &hap_cmp)
{
	auto tk = ia::trek::create_new(rov, g.get_hap_count(), is_tangled);

	const std::set<qt::up_t<qt::u32>> &reversals = hap_cmp.reversals;
	const std::set<qt::up_t<qt::u32>> &matches = hap_cmp.matches;
	const std::set<qt::up_t<qt::u32>> &mismatches = hap_cmp.mismatches;

	std::set<qt::u32> variant_refs;
	std::set<qt::u32> matches_ref;

	// always add the ref hap to the matches_ref
	matches_ref.insert(ref_h_idx);

	// auto start_t = std::chrono::steady_clock::now();

	// std::thread t1(
	//	[&matches, &matches_ref, ref_h_idx]()
	//	{
	//		for (auto [ha, hb] : matches) {
	//			if (ha == ref_h_idx)
	//				matches_ref.insert(hb);
	//			else if (hb == ref_h_idx)
	//				matches_ref.insert(ha);
	//		}
	//	});

	// std::thread t2(
	//	[&mismatches, &variant_refs, ref_h_idx]()
	//	{
	//		for (auto [ha, hb] : mismatches) {
	//			if (ha == ref_h_idx)
	//				variant_refs.insert(hb);
	//			else if (hb == ref_h_idx)
	//				variant_refs.insert(ha);
	//		}
	//	});

	for (auto [ha, hb] : matches) {
		if (ha == ref_h_idx)
			matches_ref.insert(hb);
		else if (hb == ref_h_idx)
			matches_ref.insert(ha);
	}

	for (auto [ha, hb] : mismatches) {
		if (ha == ref_h_idx)
			variant_refs.insert(hb);
		else if (hb == ref_h_idx)
			variant_refs.insert(ha);
	}

	// auto end_t = std::chrono::steady_clock::now();
	// auto elapsed = std::chrono::duration_cast<std::chrono::microseconds>(
	//	end_t - start_t);
	// INFO("Identified matches and mismatches in {} us", elapsed.count());

	// start_t = std::chrono::steady_clock::now();

	// t2.join();
	cxt_idx_t vci = find_variant_mat_contexts(mat_set, variant_refs);

	// end_t = std::chrono::steady_clock::now();
	// elapsed = std::chrono::duration_cast<std::chrono::microseconds>(
	//	end_t - start_t);
	// INFO("Identified variant contexts in {} us", elapsed.count());

	// start_t = std::chrono::steady_clock::now();

	// t1.join();
	gen_trip(g, vci, ref_h_idx, h2l, hap_itns, sorted_vertices, matches_ref,
		 mat_set.filter, tk);

	// end_t = std::chrono::steady_clock::now();
	// elapsed = std::chrono::duration_cast<std::chrono::microseconds>(
	//	end_t - start_t);

	// INFO("Generated trip in {} us", elapsed.count());

	if (tk.has_data())
		return tk;
	else
		return std::nullopt;
}

}; // namespace ita::trip
