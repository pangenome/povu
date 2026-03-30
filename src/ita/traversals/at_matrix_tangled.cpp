#include <set>
#include <vector>

#include <meza/pool/split.hpp> // for matrix_pool

#include "ita/traversals/at_matrix.hpp" // for matrix_pool
#include "ita/traversals/untangle.hpp"	// for aln_chain, chain_link
#include "ita/variation/rov.hpp"	// for RoV
#include "povu/common/core.hpp"		// for pt
#include "povu/graph/bidirected.hpp"	// for VG
#include "quilt/types.hpp"

namespace ita::at_matrix::tangled
{
namespace lq = liteseq;

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

	std::vector<qt::u8> ref_row(J, 0);

	for (qt::u32 step_idx : at) {
		const lq::ref_walk *h_w =
			g.get_ref_vec(h_idx)->walk; // the hap walk
		lq::strand lq_s = h_w->strands[step_idx];
		pt::u32 v_id = h_w->v_ids[step_idx];
		pt::u32 sorted_pos = rov->get_sorted_pos(v_id);

		qt::u8 val = (lq_s == lq::strand::STRAND_FWD) ? 1 : 2;

		ref_row[sorted_pos] = val;
	}

	// copy contents of ref row to all other rows
	for (pt::u32 i{}; i < I; i++)
		for (pt::u32 j{}; j < J; j++)
			ref_mat.set_value(i, j, ref_row[j]);
}

void populate_filter(const bd::VG &g, const ir::RoV *rov, pt::u32 h_idx,
		     const allele_traversal &at, ov_mat_t &filter_mat)
{
	if (at.empty())
		return;

	const qt::u32 J = filter_mat.cols();

	std::vector<qt::u8> ref_row(J, 0);

	for (qt::u32 step_idx : at) {
		const lq::ref_walk *h_w =
			g.get_ref_vec(h_idx)->walk; // the hap walk
		lq::strand lq_s = h_w->strands[step_idx];
		pt::u32 v_id = h_w->v_ids[step_idx];
		pt::u32 sorted_pos = rov->get_sorted_pos(v_id);

		qt::u8 val = (lq_s == lq::strand::STRAND_FWD) ? 1 : 2;

		filter_mat.set_value(h_idx, sorted_pos, val);
	}
}

void from_tangled(const bd::VG &g, const ir::RoV *rov,
		  const std::set<pt::u32> &to_call_ref_ids,
		  meza::pool::split::matrix_pool<qt::u8> &ov_pool,
		  const aln_chain &ac, rov_job_batch &batch)
{
	const std::vector<itinerary> &hap_itns = ac.hap_itns;

	pt::u32 I = g.get_hap_count();
	pt::u32 J = rov->get_vertex_count();

	rov_job j{rov, I, hap_itns};

	for (pt::u32 ref_h_idx : to_call_ref_ids) {
		if (!pv_cmp::contains(ac.all_chains, ref_h_idx))
			continue;

		ita::traversals::untangle::chain c =
			ac.all_chains.at(ref_h_idx);

		for (const auto &[ref_loop_no, links] : c.loop2ats) {

			std::map<qt::u32, qt::u32> loop_pairing;

			// -----------------------------------------------------
			// populate the ref matrix with the row of the ref hap
			// that is linked to the ref hap in the current loop
			// -----------------------------------------------------

			auto ref_mat = ov_pool.alloc_ov_matrix<std::string,
							       std::string>(
				I, J,
				meza::pool::split::pool_region::Reference);

			const allele_traversal &ref_at =
				hap_itns.at(ref_h_idx).at(ref_loop_no);

			loop_pairing.insert({ref_h_idx, ref_loop_no});

			populate_ref(g, rov, ref_h_idx, ref_at, I, J, ref_mat);

			// -----------------------------------------------------
			// filter
			//
			// -----------------------------------------------------

			auto filter_mat = ov_pool.alloc_ov_matrix<std::string,
								  std::string>(
				I, J, meza::pool::split::pool_region::Filter);

			for (qt::u32 alt_h_idx{}; alt_h_idx < I; alt_h_idx++) {

				if (alt_h_idx == ref_h_idx) {
					populate_filter(g, rov, ref_h_idx,
							ref_at, filter_mat);
					continue;
				}

				std::optional<qt::u32> opt_alt_loop_no =
					ac.get_alt_loop_no(ref_h_idx, alt_h_idx,
							   ref_loop_no);

				if (!opt_alt_loop_no)
					continue;

				loop_pairing.insert(
					{alt_h_idx, *opt_alt_loop_no});

				const allele_traversal &alt_at =
					hap_itns.at(alt_h_idx).at(
						*opt_alt_loop_no);

				populate_filter(g, rov, alt_h_idx, alt_at,
						filter_mat);
			}

			// -----------------------------------------------------
			// xor
			//
			// -----------------------------------------------------

			auto xor_mat = ov_pool.alloc_ov_matrix<std::string,
							       std::string>(
				I, J, meza::pool::split::pool_region::Xor);

			// -----------------------------------------------------
			//
			//
			// -----------------------------------------------------

			qt::u32 filter_size = filter_mat.base().size();

			mat3 m{std::move(ref_mat),
			       std::move(filter_mat),
			       std::move(xor_mat),
			       batch.get_pool_j_offset(),
			       I,
			       J,
			       loop_pairing};

			auto hap_slices =
				std::vector<std::optional<ia::hap_slice>>(I);
			tangle_info ti{ref_loop_no, std::move(hap_slices)};

			mat3_item item{std::move(m), std::move(ti)};

			j.add_item(ref_h_idx, std::move(item));

			qt::u32 offset =
				batch.get_pool_j_offset() + filter_size;
			batch.set_pool_j_offset(offset);
		}
	}

	batch.add(std::move(j));
}

} // namespace ita::at_matrix::tangled
