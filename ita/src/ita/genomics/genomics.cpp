#include "ita/genomics/genomics.hpp"

#include <algorithm> // for min, max
#include <cstddef>   // for size_t
#include <cstdlib>   // for std::max, exit, EXIT_FAILURE
#include <optional>  // for optional
#include <utility>   // for move

#include <dynamo/dynamo.hpp>  // for dynamic_interval_tree
#include <meza/pool/pool.hpp> // for matrix_pool

#include "ita/convolutions/convolutions.hpp" // for run_convs
#include "ita/genomics/allele.hpp"	     // for Exp, comp_itineraries
#include "ita/genomics/vcf.hpp"		     // for VcfRecIdx, gen_vcf_records
#include "ita/traversals/at_matrix.hpp"	     // rov_matrix_set
#include "ita/variation/rov.hpp"	     // for RoV, gen_rov
#include "ita/variation/sne.hpp"	     // for sne
#include "povu/common/app.hpp"		     // for config
#include "povu/common/core.hpp"		     // for pt, idx_t, id_t
#include "povu/common/log.hpp"		     // for ERR
#include "quilt/types.hpp"

namespace ita::genomics
{
using pool_t = meza::pool::pool<qt::u8, qt::u32>;

void find_inversions_new(const bd::VG &g,
			 const std::set<pt::id_t> &to_call_ref_ids,
			 std::map<pt::u32, std::vector<ia::inv_slice>> &res)
{
	auto foo = [](const liteseq::ref_walk *rw,
		      const std::vector<pt::u32> &positions, pt::u32 u,
		      pt::u32 v) -> pt::op_t<pt::u32>
	{
		if (rw->strands[positions[u]] == liteseq::strand::STRAND_FWD)
			return {positions[u], positions[v]};
		else
			return {positions[v], positions[u]};
	};
	using inv_it =
		dynamo::dynamic_interval_tree<pt::u32, pt::op_t<pt::u32>>;

	std::map<pt::u32, inv_it> hap_idx_to_it;
	for (pt::u32 hap_idx : to_call_ref_ids) {
		std::vector<ia::inv_slice> &inv_slices = res[hap_idx];
		dynamo::dynamic_interval_tree<pt::u32, pt::op_t<pt::u32>> t;
		const liteseq::ref_walk *rw = g.get_ref_vec(hap_idx)->walk;

		for (pt::u32 v_idx{}; v_idx < g.vtx_count(); v_idx++) {
			const std::vector<pt::u32> &positions =
				g.get_vertex_ref_idxs(v_idx, hap_idx);

			pt::u32 N = positions.size();

			if (N < 2)
				continue;

			std::optional<pt::op_t<pt::u32>> d{std::nullopt};
			for (pt::u32 i{}; i < N - 1; i++) {
				if (rw->strands[positions[i]] !=
				    rw->strands[positions[i + 1]]) {
					d = pt::op_t<pt::u32>{i, i + 1};
					break;
				}
			}

			if (d == std::nullopt)
				continue;

			auto [u, v] = *d;
			auto [fwd_pos, rev_pos] = foo(rw, positions, u, v);

			t.add({fwd_pos, fwd_pos, {fwd_pos, rev_pos}});
		}
		t.commit();

		// get covered regions
		std::vector<std::pair<pt::u32, pt::u32>> r =
			t.covered_regions();

		for (auto [s, e] : r) {
			pt::u32 len = e - s + 1;
			pt::slice fwd_slice{s, len};
			auto [_, s_rev] = t.values_at(s).front();
			pt::u32 rev_start = s_rev - len + 1;
			pt::slice_t rev_slice{rev_start, len};

			inv_slices.emplace_back(
				ia::inv_slice{rw, hap_idx, fwd_slice.start(),
					      rev_slice.start(), len});
		}
	}
}

void comp_expeditions(const bd::VG &g, std::vector<ir::RoV> &all_rovs,
		      pt::idx_t start, pt::idx_t count,
		      const std::set<pt::id_t> &to_call_ref_ids, pool_t &p,
		      ita::at_matrix::rov_job_batch &batch,
		      std::vector<ia::trek> &treks)
{
	const std::size_t N = all_rovs.size();
	for (pt::idx_t i{start}; i < start + count && i < N; i++) {
		bool last = i == ((start + count) - 1) || i == (N - 1);
		bool drain = last;

		const ir::RoV *rov = &all_rovs[i];
		ita::convolutions::populate_trips(g, rov, to_call_ref_ids, p,
						  batch, treks, drain);

		p.reset_depth_matrix(); // reset depth matrix pool for next RoV
	}

	return;
}

void gen_vcf_rec_map(const std::vector<pvst::Tree> &pvsts, bd::VG &g,
		     const std::set<pt::id_t> &to_call_ref_ids,
		     pbq::bounded_queue<iv::VcfRecIdx> &q,
		     const core::config &app_config)
{
	if (app_config.verbosity() > 0)
		INFO("Generating regions of variation (RoVs)");

	// Parse genomic region if specified
	std::optional<ir::genomic_region> region = std::nullopt;
	if (app_config.has_genomic_region()) {
		region = ir::parse_genomic_region(
			app_config.get_genomic_region().value());
		if (!region.has_value()) {
			PL_ERR("Failed to parse genomic region");
			std::exit(EXIT_FAILURE);
		}
	}

	std::vector<ir::RoV> all_rovs =
		ir::gen_rov(pvsts, g, to_call_ref_ids, region);

	const std::size_t CHUNK_SIZE = app_config.get_chunk_size();
	const std::size_t N = all_rovs.size();
	pt::u32 total_chunks = (N + CHUNK_SIZE - 1) / CHUNK_SIZE;

	if (app_config.verbosity() > 0)
		INFO("Processing chunks. Chunk count: {}", total_chunks);

	ise::pin_cushion pc;
	std::vector<ia::trek> treks;
	std::vector<ist::st> i_trees;
	std::map<pt::u32, std::vector<ia::inv_slice>> inv_slices;

	// meza::pool::split::matrix_pool<qt::u8> &ov_pool =
	//	meza::pool::split::matrix_pool<qt::u8>::init();

	auto p = meza::pool::pool<qt::u8, qt::u32>();

	// auto ov_pool =
	// meza::pool::matrix_pool<qt::u8>::create_from_megabytes();

	// ov_pool.cuda_setup_haps_xor();
	//  auto &ov_pool = meza::matrix_pool::matrix_pool<qt::u8>::init();

	// 4 bytes per u32 value
	// (1024*1024) / 4 = 262,144
	// 262144 u32 values are ~1M
	// 1024 values of u32 are 1M
	// 2,621,440 u32 values are ~10M
	// constexpr std::size_t target_bytes = 10ull * 1024 * 1024; // 10 MiB
	// constexpr std::size_t depth_matrix_pool_size =
	//	target_bytes / sizeof(qt::u32);

	// meza::pool::joint::joint_pool<qt::u32> dm_pool =
	//	meza::pool::joint::joint_pool<qt::u32>::init(
	//		depth_matrix_pool_size);

	ita::at_matrix::rov_job_batch batch;

	try {
		for (pt::idx_t base{}; base < N; base += CHUNK_SIZE) {
			pt::u32 end = std::min(base + CHUNK_SIZE, N);
			pt::u32 count = end - base;
			pt::u32 chunk_num = (base / CHUNK_SIZE) + 1;

			if (app_config.verbosity() > 0)
				INFO("\t{}/{}", chunk_num, total_chunks);

			comp_expeditions(g, all_rovs, base, count,
					 to_call_ref_ids, p, batch, treks);

			iv::VcfRecIdx rs = iv::gen_vcf_records(
				g, treks, i_trees, inv_slices);

			if (!q.push(std::move(rs)))
				break; // queue was closed early

			treks.clear();
			if (chunk_num == CHUNK_SIZE)
				i_trees.clear();
		}

		// inversions
		{
			if (app_config.verbosity() > 0)
				INFO("Processing inversions");

			find_inversions_new(g, to_call_ref_ids, inv_slices);

			iv::VcfRecIdx rs = iv::gen_vcf_records(
				g, treks, i_trees, inv_slices);
			q.push(std::move(rs));
		}
	}
	catch (...) {
		q.close(); // make sure consumers wake up on errors
		throw;
	}

	q.close(); // we're done
}
} // namespace ita::genomics
