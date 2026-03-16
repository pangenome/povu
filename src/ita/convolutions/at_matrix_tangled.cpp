#include <optional>
#include <set>
#include <vector>

#include <convo/matrix.hpp> // for ref_matrix, depth_matrix, at_matrix

#include "ita/convolutions/at_matrix.hpp" // for matrix_pool
#include "ita/genomics/allele.hpp"	  // for hap_slice
#include "ita/traversals/traversals.hpp"
#include "ita/traversals/untangle.hpp"
#include "ita/variation/rov.hpp" // for RoV
#include "liteseq/types.h"

#include "povu/common/core.hpp"	     // for pt
#include "povu/graph/bidirected.hpp" // for VG
#include "quilt/types.hpp"

namespace ita::at_matrix::tangled
{
namespace lq = liteseq;

// using ita::at_matrix::matrix_pool;
using ita::traversals::traversals::allele_traversal;
using ita::traversals::traversals::itinerary;
using ita::traversals::untangle::aln_chain;
using ita::traversals::untangle::chain_link;

using meza::matrix::depth_matrix;

void populate_ref(const bd::VG &g, const ir::RoV *rov, pt::u32 h_idx,
		  const allele_traversal &at, const pt::u32 I, const pt::u32 J,
		  ov_mat_t &ref_mat)
{
	if (at.empty())
		return;

	pt::u32 lower_bound = at.front();
	pt::u32 upper_bound = at.back();

	std::vector<qt::u8> ref_row(J, 0);

	for (pt::u32 j{}; j < J; j++) {
		pt::u32 v_id = rov->get_sorted_vertex(j);

		pt::u32 v_idx = g.v_id_to_idx(v_id);
		const std::vector<pt::idx_t> &unfiltered_step_idxs =
			g.get_vertex_ref_idxs(v_idx, h_idx);

		std::vector<pt::idx_t> in_range_step_idxs;
		for (auto step_idx : unfiltered_step_idxs)
			if (step_idx >= lower_bound && step_idx <= upper_bound)
				in_range_step_idxs.push_back(step_idx);

		pt::u32 depth = in_range_step_idxs.size();

		if (depth == 0)
			continue;

		const lq::ref_walk *h_w = g.get_ref_vec(h_idx)->walk;
		pt::u32 k = in_range_step_idxs[0]; // index in the hap walk

		lq::strand lq_s = h_w->strands[k];
		// if (lq_s == lq::strand::STRAND_FWD)
		//	ref_mat.set_value(h_idx, rov.get_sorted_pos(v_id), 1);
		// else if (lq_s == lq::strand::STRAND_REV)
		//	ref_mat.set_value(h_idx, rov.get_sorted_pos(v_id), 2);

		ref_row[rov->get_sorted_pos(v_id)] =
			(lq_s == lq::strand::STRAND_FWD) ? 1 : 2;
	}

	// copy contents of ref row to all other rows
	for (pt::u32 i{}; i < I; i++) {
		for (pt::u32 j{}; j < J; j++) {
			ref_mat.set_value(i, j, ref_row[j]);
		}
	}
}

void populate_filter(const bd::VG &g, const ir::RoV *rov, pt::u32 h_idx,
		     const allele_traversal &at, const pt::u32 I,
		     const pt::u32 J, ov_mat_t &filter_mat)
{
	if (at.empty())
		return;

	pt::u32 lower_bound = at.front();
	pt::u32 upper_bound = at.back();

	std::vector<qt::u8> ref_row(J, 0);

	for (pt::u32 j{}; j < J; j++) {
		pt::u32 v_id = rov->get_sorted_vertex(j);

		pt::u32 v_idx = g.v_id_to_idx(v_id);
		const std::vector<pt::idx_t> &unfiltered_step_idxs =
			g.get_vertex_ref_idxs(v_idx, h_idx);

		std::vector<pt::idx_t> in_range_step_idxs;
		for (auto step_idx : unfiltered_step_idxs)
			if (step_idx >= lower_bound && step_idx <= upper_bound)
				in_range_step_idxs.push_back(step_idx);

		pt::u32 depth = in_range_step_idxs.size();

		if (depth == 0)
			continue;

		const lq::ref_walk *h_w = g.get_ref_vec(h_idx)->walk;
		pt::u32 k = in_range_step_idxs[0]; // index in the hap walk

		lq::strand lq_s = h_w->strands[k];
		if (lq_s == lq::strand::STRAND_FWD)
			filter_mat.set_value(h_idx, rov->get_sorted_pos(v_id),
					     1);
		else if (lq_s == lq::strand::STRAND_REV)
			filter_mat.set_value(h_idx, rov->get_sorted_pos(v_id),
					     2);

		// ref_row[rov.get_sorted_pos(v_id)] =
		//	(lq_s == lq::strand::STRAND_FWD) ? 1 : 2;
	}

	// copy contents of ref row to all other rows
	// for (pt::u32 i{}; i < I; i++) {
	//	for (pt::u32 j{}; j < J; j++) {
	//		ref_mat.set_value(i, j, ref_row[j]);
	//	}
	// }
}

template <typename Derived, typename MatrixType, typename T, typename U,
	  typename W>
void fill_row(const bd::VG &g, const ir::RoV &rov, pt::u32 h_idx,
	      const allele_traversal &at, const pt::u32 I, const pt::u32 J,
	      meza::matrix::matrix_wrapper<Derived, MatrixType, T, U, W> &m)
{
	if (at.empty())
		return;

	pt::u32 lower_bound = at.front();
	pt::u32 upper_bound = at.back();

	for (pt::u32 j{}; j < J; j++) {
		pt::u32 v_id = rov.get_sorted_vertex(j);

		pt::u32 v_idx = g.v_id_to_idx(v_id);
		const std::vector<pt::idx_t> &unfiltered_step_idxs =
			g.get_vertex_ref_idxs(v_idx, h_idx);

		std::vector<pt::idx_t> in_range_step_idxs;
		for (auto step_idx : unfiltered_step_idxs)
			if (step_idx >= lower_bound && step_idx <= upper_bound)
				in_range_step_idxs.push_back(step_idx);

		pt::u32 depth = in_range_step_idxs.size();

		if (depth == 0)
			continue;

		const lq::ref_walk *h_w = g.get_ref_vec(h_idx)->walk;
		pt::u32 k = in_range_step_idxs[0]; // index in the hap walk

		lq::strand lq_s = h_w->strands[k];
		if (lq_s == lq::strand::STRAND_FWD)
			m.set_value(h_idx, rov.get_sorted_pos(v_id), 1);
		else if (lq_s == lq::strand::STRAND_REV)
			m.set_value(h_idx, rov.get_sorted_pos(v_id), 2);
	}
}

void from_tangled(const bd::VG &g, const ir::RoV *rov,
		  const std::set<pt::u32> &to_call_ref_ids,
		  meza::matrix_pool::matrix_pool<qt::u8> &ov_pool,
		  const aln_chain &ac, rov_job_batch &batch)
{

	const std::vector<itinerary> &hap_itns = ac.hap_itns;

	std::vector<mat3_item> items;

	pt::u32 I = g.get_hap_count();
	pt::u32 J = rov->get_vertex_count();

	qt::u32 j_off = batch.pool_j_offset;

	rov_job j{rov, hap_itns};

	for (pt::u32 ref_h_idx : to_call_ref_ids) {
		if (!pv_cmp::contains(ac.all_chains, ref_h_idx))
			continue;

		ita::traversals::untangle::chain c =
			ac.all_chains.at(ref_h_idx);

		for (const auto &[ref_loop_no, links] : c.loop2ats) {

			// -----------------------------------------------------
			// populate the ref matrix with the row of the ref hap
			// that is linked to the ref hap in the current loop
			// -----------------------------------------------------

			auto ref_mat = ov_pool.alloc_ov_matrix<std::string,
							       std::string>(
				I, J,
				meza::matrix_pool::pool_region::Reference);

			const allele_traversal &ref_at =
				hap_itns.at(ref_h_idx).at(ref_loop_no);

			populate_ref(g, rov, ref_h_idx, ref_at, I, J, ref_mat);

			// -----------------------------------------------------
			// filter
			//
			// -----------------------------------------------------

			auto filter_mat = ov_pool.alloc_ov_matrix<std::string,
								  std::string>(
				I, J, meza::matrix_pool::pool_region::Filter);

			auto hap_slices =
				std::vector<std::optional<ia::hap_slice>>(I);

			for (const ita::traversals::untangle::chain_link &link :
			     links) {
				auto [_, __, ___, alt_h_idx, alt_loop_no] =
					link;

				const allele_traversal &alt_at =
					hap_itns.at(alt_h_idx).at(alt_loop_no);

				populate_filter(g, rov, alt_h_idx, alt_at, I, J,
						filter_mat);
			}

			// -----------------------------------------------------
			// xor
			//
			// -----------------------------------------------------

			auto xor_mat = ov_pool.alloc_ov_matrix<std::string,
							       std::string>(
				I, J, meza::matrix_pool::pool_region::Xor);

			// -----------------------------------------------------
			//
			//
			// -----------------------------------------------------

			mat3 m{std::move(ref_mat),
			       std::move(filter_mat),
			       std::move(xor_mat),
			       j_off,
			       I,
			       J};

			tangle_info ti{ref_loop_no, std::move(hap_slices)};

			mat3_item item{std::move(m), std::move(ti)};

			j.add_item(ref_h_idx, std::move(item));

			// items.emplace_back(std::move(item));

			// rov_job j{
			//	rov,
			//	batch.pool_j_offset,
			//	ref_h_idx,
			//	std::move(item),
			// };

			// items.clear();

			j_off += J;
		}
	}

	batch.add(std::move(j));
	batch.pool_j_offset = j_off;
}

// rov_matrix_pool
// init_tangled_depth_matrices(const bd::VG &g, const ir::RoV &rov,
//			    const std::set<pt::u32> &to_call_ref_ids,
//			    const aln_chain &ac)
// {
//	rov_matrix_pool res{rov, true, ac.hap_itns};

//	pt::u32 I = g.get_hap_count();
//	pt::u32 J = rov.get_vertex_count(); // cols
//	const std::vector<pt::u32> &sorted_vertices = rov.get_sorted_vertices();
//	const std::vector<itinerary> &hap_itns = ac.hap_itns;

//	for (pt::u32 ref_h_idx : to_call_ref_ids) {
//		if (!pv_cmp::contains(ac.all_chains, ref_h_idx))
//			continue;

//		ita::traversals::untangle::chain c =
//			ac.all_chains.at(ref_h_idx);

//		for (const auto &[ref_loop_no, links] : c.loop2ats) {

//			// -----------------------------------------------------
//			// populate the ref matrix with the row of the ref hap
//			// that is linked to the ref hap in the current loop
//			// -----------------------------------------------------
//			auto ref_matrix = meza::matrix::ref_matrix{I, J};
//			auto cpy = sorted_vertices;
//			ref_matrix.add_col_names(std::move(cpy));

//			const allele_traversal &ref_at =
//				hap_itns.at(ref_h_idx).at(ref_loop_no);

//			fill_row(g, rov, ref_h_idx, ref_at, I, J, ref_matrix);

//			// -----------------------------------------------------
//			//
//			//
//			// -----------------------------------------------------

//			auto filter_matrix = meza::matrix::at_matrix{I, J};
//			cpy = sorted_vertices;
//			filter_matrix.add_col_names(std::move(cpy));

//			auto hap_slices =
//				std::vector<std::optional<ia::hap_slice>>(I);

//			// duplicate the filter matrix
//			auto result_matrix = filter_matrix;

//			fill_row(g, rov, ref_h_idx, ref_at, I, J,
//				 filter_matrix);

//			for (const ita::traversals::untangle::chain_link &link :
//			     links) {
//				auto [_, __, ___, alt_h_idx, alt_loop_no] =
//					link;

//				const allele_traversal &alt_at =
//					hap_itns.at(alt_h_idx).at(alt_loop_no);

//				fill_row(g, rov, alt_h_idx, alt_at, I, J,
//					 filter_matrix);
//			}

//			std::map<pt::u32, meza::matrix::ref_matrix>
//				ref_matrices{
//					{ref_h_idx, std::move(ref_matrix)}};

//			matrix_pool mp{std::move(ref_matrices),
//				       std::move(filter_matrix),
//				       std::move(result_matrix),
//				       std::move(hap_slices),
//				       ref_loop_no,
//				       true,
//				       I,
//				       J};

//			res.pools.emplace_back(std::move(mp));
//		}
//	}

//	return res;
// }

} // namespace ita::at_matrix::tangled
