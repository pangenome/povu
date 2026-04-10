#include "ita/convolutions/convolutions.hpp"

#include <cassert>
#include <optional>
#include <set>
#include <vector>

#include <liteseq/refs.h> // for ref_walk, ref

#include <meza/pool/pool.hpp> // for pool

#include "ita/convolutions/trip.hpp"	   // for gen_trip
#include "ita/genomics/allele.hpp"	   // for hap_slice, trek
#include "ita/traversals/at_matrix.hpp"	   // for matrix_pool, rov_matrix_set
#include "ita/traversals/depth_matrix.hpp" // for comp_depth_matrix
#include "ita/variation/rov.hpp"	   // for RoV
#include "povu/common/core.hpp"		   // for pt
#include "povu/graph/bidirected.hpp"	   // for VG
#include "quilt/types.hpp"		   // for u8, u32

namespace ita::convolutions
{
namespace lq = liteseq;
using meza::pool::hap_comp::haps_comp_set;

using pool_t = meza::pool::pool<qt::u8, qt::u32>;

void process_mat3(
	const bd::VG &g, const ir::RoV *rov, qt::u32 ref_h_idx,
	const std::vector<ita::traversals::traversals::itinerary> &hap_itns,
	const std::vector<ita::at_matrix::mat3_item> &ji, pool_t &p,
	std::vector<ia::trek> &treks)
{
	const std::vector<pt::u32> &sorted_vertices =
		rov->get_sorted_vertices();

	for (const ita::at_matrix::mat3_item &item : ji) {
		const ita::at_matrix::mat3 &mat_set = item.mats;
		const meza::pool::ov_mat_t &filter_mat = mat_set.filter;

		if (filter_mat.base().is_row_blank(ref_h_idx))
			continue;

		qt::u32 pool_offset = mat_set.j_offset;

		bool is_tangled = item.is_tangled();

		const ita::at_matrix::hap2loop &h2l = mat_set.h2l;

		const haps_comp_set &hap_cmp =
			p.handle_set(filter_mat, pool_offset);

		std::optional<ia::trek> opt_tk = ita::trip::gen_trip(
			g, rov, is_tangled, ref_h_idx, h2l, hap_itns, mat_set,
			sorted_vertices, hap_cmp);

		if (!opt_tk.has_value())
			continue;

		treks.emplace_back(std::move(*opt_tk));
	}
}

void process_batches(const bd::VG &g, const std::set<pt::u32> &to_call_ref_ids,
		     pool_t &p, ita::at_matrix::rov_job_batch &batch,
		     std::vector<ia::trek> &treks)
{
	p.run_convolutions(batch.get_pool_j_offset());

	for (const ita::at_matrix::rov_job &job : batch.get_jobs()) {
		const ir::RoV *rov = job.get_rov();

		const std::vector<ita::traversals::traversals::itinerary>
			&hap_itns = job.get_hap_itns();

		for (qt::u32 ref_h_idx : to_call_ref_ids) {
			const std::vector<ita::at_matrix::mat3_item> &job_item =
				job.get_items_for_ref2(ref_h_idx);

			process_mat3(g, rov, ref_h_idx, hap_itns, job_item, p,
				     treks);
		}
	}

	p.clear_split_pool();
	batch.reset();
}

void populate_trips(const bd::VG &g, const ir::RoV *rov,
		    const std::set<pt::u32> &to_call_ref_ids, pool_t &p,
		    ita::at_matrix::rov_job_batch &batch,
		    std::vector<ia::trek> &treks, bool drain)
{
	if (p.is_full() || drain)
		process_batches(g, to_call_ref_ids, p, batch, treks);

	ita::depth_matrix::depth_matrix dm =
		ita::depth_matrix::comp_depth_matrix(g, rov, p);

	// loads data from the depth matrix into
	ita::at_matrix::init_pool(g, rov, to_call_ref_ids, dm, p, batch);
}

} // namespace ita::convolutions
