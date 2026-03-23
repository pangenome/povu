#include "ita/convolutions/convolutions.hpp"

// #include <algorithm>
#include <cassert>
#include <optional>
#include <set>
// #include <thread>
#include <vector>

#include <meza/pool/joint.hpp> // for joint_pool
#include <meza/pool/split.hpp> // for matrix_pool

// #include <convo/pool_joint.hpp> // for joint_pool
// #include <convo/pool_split.hpp> // for matrix_pool, ov_mat_t

// #include <convo/pool.hpp> // for matrix_pool,
#include <liteseq/refs.h> // for ref_walk, ref

#include "ita/convolutions/at_matrix.hpp"    // for matrix_pool, rov_matrix_set
#include "ita/convolutions/depth_matrix.hpp" // for comp_depth_matrix
#include "ita/convolutions/trip.hpp"	     // for gen_trip
#include "ita/genomics/allele.hpp"	     // for hap_slice, trek
#include "ita/variation/rov.hpp"	     // for RoV
// #include "povu/common/constants.hpp"
#include "povu/common/core.hpp"	     // for pt
#include "povu/graph/bidirected.hpp" // for VG
#include "quilt/types.hpp"

namespace ita::convolutions
{
namespace lq = liteseq;

void comp_conv(meza::pool::split::matrix_pool<qt::u8> &ov_pool, qt::u32 I,
	       qt::u32 pool_j_offset)
{
	// copy the pool to the device before launching the kernel
	ov_pool.copy_to_device();

	ov_pool.xor_on_device(pool_j_offset);

	// copy the pool back to the host after the kernel has finished
	ov_pool.copy_to_host_thirds(meza::pool::split::pool_region::Xor);

	// ov_pool.sync_device();
}

void process_batches(const bd::VG &g, const std::set<pt::u32> &to_call_ref_ids,
		     meza::pool::split::matrix_pool<qt::u8> &ov_pool,
		     ita::at_matrix::rov_job_batch &batch,
		     std::vector<ia::trek> &treks)
{
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

			for (const ita::at_matrix::mat3_item &item : *ji) {

				const ita::at_matrix::mat3 &mat_set = item.mats;
				const meza::pool::split::ov_mat_t &filter_mat =
					mat_set.filter;
				qt::u32 pool_offset = mat_set.j_offset;

				bool is_tangled = item.is_tangled();

				const ita::at_matrix::hap2loop &h2l =
					mat_set.h2l;

				// hap2loop h2l(mat_set.loop_pairing);

				const meza::pool::split::haps_comp_set &
					hap_cmp = meza::pool::split::handle_set(
						ov_pool, filter_mat,
						pool_offset);

				std::optional<ia::trek> opt_tk =
					ita::trip::gen_trip(
						g, rov, is_tangled, ref_h_idx,
						h2l, hap_itns, mat_set,
						sorted_vertices, hap_cmp);

				if (!opt_tk.has_value())
					continue;

				treks.emplace_back(std::move(*opt_tk));
			}
		}
	}

	ov_pool.clear();
	batch.reset();
}

void populate_trips(const bd::VG &g, const ir::RoV *rov,
		    const std::set<pt::u32> &to_call_ref_ids,
		    meza::pool::joint::joint_pool<qt::u32> &dm_pool,
		    meza::pool::split::matrix_pool<qt::u8> &ov_pool,
		    ita::at_matrix::rov_job_batch &batch,
		    std::vector<ia::trek> &treks, bool drain)
{

	if (ov_pool.is_full() || drain) {
		std::cerr << "draining at " << rov->as_str() << "\n";
		process_batches(g, to_call_ref_ids, ov_pool, batch, treks);
	};

	ita::depth_matrix::depth_matrix dm =
		ita::depth_matrix::comp_depth_matrix(g, rov, dm_pool);

	// loads data from the depth matrix into
	ita::at_matrix::init_pool(g, rov, to_call_ref_ids, dm, ov_pool, batch);

	// std::cerr << "pool is full? " << ov_pool.is_full() << " j off "
	//	  << batch.pool_j_offset << "\n";
	// auto [ref_used, filter_used, xor_used] = ov_pool.used();
	// std::cerr << ref_used << " " << filter_used << " " << xor_used <<
	// "\n";

	// std::exit(1);
}

} // namespace ita::convolutions
