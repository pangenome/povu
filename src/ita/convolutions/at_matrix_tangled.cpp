#include <optional>
#include <set>
#include <vector>

#include "ita/convolutions/at_matrix.hpp" // for matrix_pool
#include "ita/genomics/allele.hpp"	  // for hap_slice
#include "ita/traversals/traversals.hpp"
#include "ita/traversals/untangle.hpp"
#include "ita/variation/rov.hpp" // for RoV
#include "liteseq/types.h"
#include "meza/matrix/matrix.hpp"
#include "povu/common/core.hpp"	     // for pt
#include "povu/graph/bidirected.hpp" // for VG

namespace ita::at_matrix::tangled
{
namespace lq = liteseq;

using ita::at_matrix::matrix_pool;
using ita::traversals::traversals::allele_traversal;
using ita::traversals::traversals::itinerary;
using ita::traversals::untangle::aln_chain;
using ita::traversals::untangle::chain_link;

using meza::matrix::depth_matrix;

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

rov_matrix_pool
init_tangled_depth_matrices(const bd::VG &g, const ir::RoV &rov,
			    const std::set<pt::u32> &to_call_ref_ids,
			    const aln_chain &ac)
{
	rov_matrix_pool res{rov, true, ac.hap_itns};

	pt::u32 I = g.get_hap_count();
	pt::u32 J = rov.get_vertex_count(); // cols
	const std::vector<pt::u32> &sorted_vertices = rov.get_sorted_vertices();
	const std::vector<itinerary> &hap_itns = ac.hap_itns;

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
			auto ref_matrix = meza::matrix::ref_matrix{I, J};
			auto cpy = sorted_vertices;
			ref_matrix.add_col_names(std::move(cpy));

			const allele_traversal &ref_at =
				hap_itns.at(ref_h_idx).at(ref_loop_no);

			fill_row(g, rov, ref_h_idx, ref_at, I, J, ref_matrix);

			// -----------------------------------------------------
			//
			//
			// -----------------------------------------------------

			auto filter_matrix = meza::matrix::at_matrix{I, J};
			cpy = sorted_vertices;
			filter_matrix.add_col_names(std::move(cpy));

			auto hap_slices =
				std::vector<std::optional<ia::hap_slice>>(I);

			// duplicate the filter matrix
			auto result_matrix = filter_matrix;

			fill_row(g, rov, ref_h_idx, ref_at, I, J,
				 filter_matrix);

			for (const ita::traversals::untangle::chain_link &link :
			     links) {
				auto [_, __, ___, alt_h_idx, alt_loop_no] =
					link;

				const allele_traversal &alt_at =
					hap_itns.at(alt_h_idx).at(alt_loop_no);

				fill_row(g, rov, alt_h_idx, alt_at, I, J,
					 filter_matrix);
			}

			std::map<pt::u32, meza::matrix::ref_matrix>
				ref_matrices{
					{ref_h_idx, std::move(ref_matrix)}};

			matrix_pool mp{std::move(ref_matrices),
				       std::move(filter_matrix),
				       std::move(result_matrix),
				       std::move(hap_slices),
				       ref_loop_no,
				       true,
				       I,
				       J};

			res.pools.emplace_back(std::move(mp));
		}
	}

	return res;
}

} // namespace ita::at_matrix::tangled
