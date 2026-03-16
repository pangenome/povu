#include <cassert>
#include <optional>

#include <convo/pool.hpp> // for matrix_pool
#include <liteseq/refs.h> // for ref_walk, ref

#include "convo/matrix.hpp"
#include "ita/convolutions/depth_matrix.hpp" // for comp_depth_matrix

#include "ita/convolutions/at_matrix.hpp" // for matrix_pool, rov_matrix_set
#include "ita/convolutions/trip.hpp"	  // for gen_trips
#include "ita/genomics/allele.hpp"	  // for hap_slice, trek
#include "ita/variation/rov.hpp"	  // for RoV
#include "povu/common/core.hpp"		  // for pt
#include "povu/common/utils.hpp"
#include "povu/graph/bidirected.hpp" // for VG
#include "quilt/types.hpp"

namespace ita::convolutions
{
namespace lq = liteseq;
using u32 = qt::u32;
using ov_mat_t = meza::matrix_view::ov_matrix<qt::u8, std::string, std::string>;

struct comp_res {
	std::map<pt::op_t<pt::u32>, std::set<pt::u32>> contexts2haps;
};

void populate_filter(const qt::u32 I, const qt::u32 J,
		     const meza::matrix::depth_matrix &dm_filter_mat,
		     ov_mat_t &filter_mat)
{
	for (qt::u32 i{}; i < I; i++) {
		for (qt::u32 j{}; j < J; j++) {
			qt::u8 val = dm_filter_mat.base().at(i, j);
			filter_mat.set_value(i, j, val);
		}
	}
}

void populate_ref(const qt::u32 I, const qt::u32 J, qt::u32 ref_i,
		  const meza::matrix::depth_matrix &dm_filter_mat,
		  ov_mat_t &ref_mat)
{
	for (qt::u32 i{}; i < I; i++) {
		for (qt::u32 j{}; j < J; j++) {
			qt::u8 val = dm_filter_mat.base().at(ref_i, j);
			ref_mat.set_value(i, j, val);
		}
	}
}

void comp_conv(meza::matrix_pool::matrix_pool<qt::u8> &ov_pool, qt::u32 I,
	       qt::u32 total_J)
{
	// copy the pool to the device before launching the kernel
	ov_pool.copy_to_device();

	ov_pool.xor_on_device(I, total_J);

	// copy the pool back to the host after the kernel has finished
	ov_pool.copy_to_host_thirds(meza::matrix_pool::pool_region::Xor);
}

std::vector<qt::op_t<qt::u32>> find_context(const ov_mat_t &ref_matrix,
					    const ov_mat_t &xor_mat, qt::u32 i)
{
	qt::u32 J = ref_matrix.cols();
	if (J == 1)
		return {};

	std::vector<qt::u32> sum_row(J, 0);
	sum_row[0] = xor_mat.base().at(i, 0);
	for (u32 j{1}; j < J; j++) {
		qt::u8 val = xor_mat.base().at(i, j);
		sum_row[j] = sum_row[j - 1] + val;
	}

	std::vector<qt::op_t<qt::u32>> bounds;
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

std::map<pt::u32, std::vector<pt::op_t<pt::u32>>>
find_context2(const ita::at_matrix::mat3 &mat_set,
	      const std::set<pt::u32> &variant_refs)
{
	const ov_mat_t &ref_mat = mat_set.ref;
	const ov_mat_t &xor_mat = mat_set.xor_result;

	std::map<pt::u32, std::vector<pt::op_t<pt::u32>>> hap_contexts;
	for (pt::u32 h_idx : variant_refs) {
		std::vector<pt::op_t<pt::u32>> contexts =
			find_context(ref_mat, xor_mat, h_idx);
		hap_contexts.emplace(h_idx, std::move(contexts));
	}

	return hap_contexts;
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

comp_res foo(const ita::at_matrix::mat3 &mat_set,
	     const std::set<pt::u32> &variant_refs)
{
	std::map<pt::u32, std::vector<pt::op_t<pt::u32>>> hap_contexts =
		find_context2(mat_set, variant_refs);

	std::map<pt::op_t<pt::u32>, std::set<pt::u32>> contexts2haps =
		contexts_to_hap_idxs(hap_contexts);

	return {contexts2haps};
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
	pt::u32 u_v_idx = g.v_id_to_idx(sorted_vertices[u]);
	pt::u32 v_v_idx = g.v_id_to_idx(sorted_vertices[v]);

	pt::u32 lower_bound = at.front();
	pt::u32 upper_bound = at.back();

	pt::u32 u_step_idx;
	for (qt::u8 step_idx : g.get_vertex_ref_idxs(u_v_idx, h_idx))
		if (step_idx >= lower_bound && step_idx <= upper_bound)
			u_step_idx = step_idx;

	pt::u32 v_step_idx;
	for (qt::u8 step_idx : g.get_vertex_ref_idxs(v_v_idx, h_idx))
		if (step_idx >= lower_bound && step_idx <= upper_bound)
			v_step_idx = step_idx;

	ptg::id_or_t u_id_or = step_idx2id_or(u_step_idx);
	ptg::id_or_t v_id_or = step_idx2id_or(v_step_idx);

	return ia::rov_boundaries{u_id_or, v_id_or};
}

pt::op_t<pt::u32>
matrix_context_to_graph_context(const std::vector<pt::u32> &sorted_vertices,
				pt::op_t<pt::u32> context)
{
	auto [u, v] = context;
	pt::u32 u_v_id = sorted_vertices[u];
	pt::u32 v_v_id = sorted_vertices[v];

	return {u_v_id, v_v_id};
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
gen_hap_slices(const bd::VG &g, const std::vector<pt::u32> &sorted_vertices,
	       pt::op_t<pt::u32> context, pt::u32 ref_h_idx,
	       const std::set<pt::u32> &h_idxs)
{
	pt::op_t<pt::u32> graph_context =
		matrix_context_to_graph_context(sorted_vertices, context);

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

void gen_trip(
	const bd::VG &g, const comp_res &cr, qt::u32 ref_h_idx,
	qt::u32 ref_loop_no,
	const std::vector<ita::traversals::traversals::itinerary> &hap_itns,
	bool is_tangled, const std::vector<pt::u32> &sorted_vertices,
	const std::set<qt::u32> &matches_ref, const ov_mat_t &filter_mat,
	ia::trek &tk)
{
	const qt::u32 I = g.get_hap_count();
	const auto &[context_to_haps] = cr;

	for (const auto &[context, haps] : context_to_haps) {

		std::optional<ia::rov_boundaries> opt_ia_cxt;
		if (is_tangled) {
			const ita::traversals::traversals::itinerary
				&ref_hap_itn = hap_itns.at(ref_h_idx);
			const ita::traversals::traversals::allele_traversal
				&at = ref_hap_itn.at(ref_loop_no);
			opt_ia_cxt = indexes_to_rov_boundaries_tangled(
				g, context, sorted_vertices, at, ref_h_idx);
		}
		else {
			opt_ia_cxt = indexes_to_rov_boundaries_no_tangle(
				g, context, sorted_vertices, ref_h_idx);
		}

		ia::rov_boundaries ia_cxt = opt_ia_cxt.value();

		auto [ref_slice, alt_set, alt_haps] = gen_hap_slices(
			g, sorted_vertices, context, ref_h_idx, haps);

		ia::cxt_to_min_rov_map d;
		d.emplace(ia_cxt,
			  ia::minimal_rov(ia_cxt, ref_slice, matches_ref,
					  alt_haps, std::move(alt_set)));

		tk.set_min_rov(ref_h_idx, std::move(d));

		for (auto h_idx : matches_ref)
			tk.add_match_ref(ref_h_idx, h_idx);

		for (pt::u32 h_idx{}; h_idx < I; h_idx++)
			if (filter_mat.base().is_row_blank(h_idx) &&
			    h_idx != ref_h_idx)
				tk.add_no_cov(ref_h_idx, h_idx);
	}
}

void process_batches(const bd::VG &g, const std::set<pt::u32> &to_call_ref_ids,
		     meza::matrix_pool::matrix_pool<qt::u8> &ov_pool,
		     ita::at_matrix::rov_job_batch &batch,
		     std::vector<ia::trek> &treks)
{
	// qt::u32 I = g.get_hap_count();
	// qt::u32 pool_j_offset = batch.pool_j_offset;
	comp_conv(ov_pool, g.get_hap_count(), batch.pool_j_offset);

	for (const ita::at_matrix::rov_job &job : batch.items) {
		const ir::RoV *rov = job.rov;
		const std::vector<pt::u32> &sorted_vertices =
			rov->get_sorted_vertices();

		const std::vector<ita::traversals::traversals::itinerary>
			&hap_itns = job.hap_itns;

		for (qt::u32 ref_h_idx : to_call_ref_ids) {
			const std::vector<ita::at_matrix::mat3_item> *ji =
				job.get_items_for_ref(ref_h_idx);

			if (ji == nullptr)
				continue;

			qt::u32 N = ji->size();
			bool is_tangled = N > 1 ? true : false;

			std::cerr << rov->as_str() << "\n";

			for (qt::u32 ref_loop_no{}; ref_loop_no < N;
			     ref_loop_no++) {
				const ita::at_matrix::mat3_item &item =
					ji->at(ref_loop_no);

				const ita::at_matrix::mat3 &mat_set = item.mats;
				const ov_mat_t &filter_mat = mat_set.filter;
				qt::u32 pool_offset = mat_set.j_offset;

				// const auto &[_, filter_mat, ___, pool_offset]
				// =	mat_set;

				std::cerr << __func__ << " pool_offset "
					  << pool_offset << "\n";

				const auto &[reversals, matches, mismatches] =
					meza::matrix_pool::handle_set(
						ov_pool, filter_mat,
						pool_offset);

				// const auto &[matches, mismatches] =
				//	ita::trip::comp_at_comparison(ov_pool,
				//				      mat_set);

				auto tk = ia::trek::create_new(
					rov, g.get_hap_count(), false);

				std::set<qt::u32> variant_refs;
				std::set<qt::u32> matches_ref;

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

				if (rov->as_str() == ">1546>1551") {
					mat_set.dbg_print();
				}

				comp_res cr = foo(mat_set, variant_refs);
				auto [contexts2haps] = cr;
				std::cerr << "contexts2haps size: "
					  << contexts2haps.size() << "\n";

				for (const auto &[context, haps] :
				     contexts2haps) {
					std::cerr << "(" << context.first
						  << ", " << context.second
						  << "): {";
					for (qt::u32 h_idx : haps)
						std::cerr << h_idx << ", ";
					std::cerr << "}\n";
				}

				gen_trip(g, cr, ref_h_idx, ref_loop_no,
					 hap_itns, is_tangled, sorted_vertices,
					 matches_ref, mat_set.filter, tk);

				if (tk.has_data()) {
					std::cerr << "trek for ref "
						  << ref_h_idx
						  << " has data, "
						     "adding to "
						     "treks\n";

					treks.emplace_back(std::move(tk));
				}
				else {
					std::cerr << ref_loop_no
						  << " trek is blank\n";

					std::cerr << "variant refs: \n";
					std::cerr << pu::concat_with(
							     variant_refs, ',')
						  << "\n";

					mat_set.dbg_print();
					std::cerr << rov->as_str()
						  << " early exit\n";
					std::exit(1);
				}
			}
		}
	}

	ov_pool.reset();
	batch.reset();
}

// void run_drain(const bd::VG &g, const ir::RoV *rov, const qt::u32 I,
//	       const std::set<pt::u32> &to_call_ref_ids,
//	       meza::matrix_pool::matrix_pool<qt::u8> &ov_pool,
//	       ita::at_matrix::rov_job_batch &batch,
//	       std::vector<ia::trek> &treks)
// {
//	qt::u32 pool_j_offset = batch.pool_j_offset;
//	comp_conv(ov_pool, I, pool_j_offset);

//	for (const ita::at_matrix::rov_job &job : batch.items) {

//		const std::vector<pt::u32> &sorted_vertices =
//			job.rov->get_sorted_vertices();

//		for (qt::u32 ref_h_idx : to_call_ref_ids) {
//			const std::vector<ita::at_matrix::mat3_item> *ji =
//				job.get_items_for_ref(ref_h_idx);

//			if (ji == nullptr)
//				continue;

//			for (const ita::at_matrix::mat3_item &item : *ji) {
//				const ita::at_matrix::mat3 &mat_set = item.mats;
//				const auto &[matches, mismatches] =
//					ita::trip::comp_at_comparison(ov_pool,
//								      mat_set);

//				auto tk = ia::trek::create_new(
//					rov, g.get_hap_count(), false);

//				for (qt::u32 ref_h_idx : to_call_ref_ids) {
//					std::set<qt::u32> variant_refs;
//					std::set<qt::u32> matches_ref;

//					for (auto [ha, hb] : matches) {
//						if (ha == ref_h_idx)
//							matches_ref.insert(hb);
//						else if (hb == ref_h_idx)
//							matches_ref.insert(ha);
//					}

//					for (auto [ha, hb] : mismatches) {
//						if (ha == ref_h_idx)
//							variant_refs.insert(hb);
//						else if (hb == ref_h_idx)
//							variant_refs.insert(ha);
//					}

//					comp_res cr =
//						foo(mat_set, variant_refs);

//					gen_trip(g, I, cr, ref_h_idx,
//						 sorted_vertices, matches_ref,
//						 mat_set.filter, tk);

//					if (tk.has_data()) {
//						std::cerr << "trek for ref "
//							  << ref_h_idx
//							  << " has data, "
//							     "adding to "
//							     "treks\n";

//						treks.emplace_back(
//							std::move(tk));
//					}
//				}
//			}
//		}
//	}

//	ov_pool.reset();
//	batch.reset();
// }

void populate_trips(const bd::VG &g, const ir::RoV *rov,
		    const std::set<pt::u32> &to_call_ref_ids,
		    meza::matrix_pool::joint_pool<qt::u32> &dm_pool,
		    meza::matrix_pool::matrix_pool<qt::u8> &ov_pool,
		    ita::at_matrix::rov_job_batch &batch,
		    std::vector<ia::trek> &treks, bool drain)
{
	if (drain) {
		std::cerr << "draining at " << rov->as_str() << "\n";
		process_batches(g, to_call_ref_ids, ov_pool, batch, treks);
	};

	ita::depth_matrix::depth_matrix dm =
		ita::depth_matrix::comp_depth_matrix(g, rov, dm_pool);

	// loads data from the depth matrix into
	ita::at_matrix::init_pool(g, rov, to_call_ref_ids, dm, ov_pool, batch);
}

// void comp_expedition(const bd::VG &g, const ir::RoV *rov, bool last,
//		     const std::set<pt::u32> &to_call_ref_ids,
//		     meza::matrix_pool::matrix_pool<qt::u8> &ov_pool,
//		     ita::at_matrix::all_mat_sets &mat_sets,
//		     std::vector<ia::trek> &treks)
// {
//	using ita::at_matrix::matrix_pool;
//	using ita::at_matrix::rov_matrix_pool;

//	rov_matrix_pool rov_mp =
//		ita::at_matrix::init_depth_matrices(g, rov, to_call_ref_ids);

//	if (rov_mp.tangled)
//		return;

//	qt::u32 I = rov_mp.pools.front().I;
//	qt::u32 J = rov_mp.pools.front().J;

//	bool is_pool_full = !ov_pool.can_allocate(
//		I, J * 3, meza::matrix_view::layout::DenseRowMajor);

//	if (is_pool_full || last) {

//		std::cerr << __func__ << " emptying pool " << "\n";

//		comp_conv(ov_pool, I, mat_sets.pool_j_offset);

//		for (const auto &mat_set : mat_sets.mat_sets) {
//			const auto &[matches, mismatches] =
//				ita::trip::comp_at_comparison(ov_pool, mat_set);

//			auto tk = ia::trek::create_new(rov, g.get_hap_count(),
//						       false);

//			const std::vector<pt::u32> &sorted_vertices =
//				mat_set.rov->get_sorted_vertices();

//			for (qt::u32 ref_h_idx : to_call_ref_ids) {
//				std::set<qt::u32> variant_refs;
//				std::set<qt::u32> matches_ref;

//				for (auto [ha, hb] : matches) {
//					if (ha == ref_h_idx)
//						matches_ref.insert(hb);
//					else if (hb == ref_h_idx)
//						matches_ref.insert(ha);
//				}

//				for (auto [ha, hb] : mismatches) {
//					if (ha == ref_h_idx)
//						variant_refs.insert(hb);
//					else if (hb == ref_h_idx)
//						variant_refs.insert(ha);
//				}

//				comp_res cr = foo(mat_set, variant_refs);

//				gen_trip(g, I, cr, ref_h_idx, sorted_vertices,
//					 matches_ref, mat_set.filter, tk);

//				if (tk.has_data()) {
//					std::cerr << "trek for ref "
//						  << ref_h_idx
//						  << " has data, adding to "
//						     "treks\n";

//					treks.emplace_back(std::move(tk));
//				}
//			}
//		}

//		ov_pool.reset();
//		mat_sets.reset();
//	}

//	for (const matrix_pool &pool : rov_mp.pools) {

//		const meza::matrix::depth_matrix &dm_filter_mat =
//			pool.filter_matrix;

//		for (qt::u32 ref_h_idx : to_call_ref_ids) {

//			auto ref_mat = ov_pool.alloc_ov_matrix<std::string,
//							       std::string>(
//				I, J,
//				meza::matrix_pool::pool_region::Reference);

//			auto filter_mat = ov_pool.alloc_ov_matrix<std::string,
//								  std::string>(
//				I, J, meza::matrix_pool::pool_region::Filter);

//			auto xor_mat = ov_pool.alloc_ov_matrix<std::string,
//							       std::string>(
//				I, J, meza::matrix_pool::pool_region::Xor);

//			populate_filter(I, J, dm_filter_mat, filter_mat);

//			populate_ref(I, J, ref_h_idx, dm_filter_mat, ref_mat);

//			mat_sets.add({ref_mat, filter_mat, xor_mat,
//				      mat_sets.pool_j_offset, rov});

//			mat_sets.pool_j_offset += J;
//		}
//	}

//	return;
// }

} // namespace ita::convolutions
