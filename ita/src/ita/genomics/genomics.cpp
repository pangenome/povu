#include "ita/genomics/genomics.hpp"

#include <algorithm> // for min, max
#include <cstddef>   // for size_t
#include <cstdlib>   // for std::max, exit, EXIT_FAILURE
#include <optional>  // for optional
#include <utility>   // for move

#include <dynamo/dynamo.hpp> // for dynamic_interval_tree
#include <log/log.h>
#include <meza/pool/pool.hpp> // for matrix_pool
#include <quilt/app.hpp>      // for config
#include <quilt/types.hpp>    // for qt

#include "ita/convolutions/convolutions.hpp" // for run_convs
#include "ita/genomics/allele.hpp"	     // for Exp, comp_itineraries
#include "ita/genomics/vcf.hpp"		     // for VcfRecIdx, gen_vcf_records
#include "ita/traversals/at_matrix.hpp"	     // rov_matrix_set
#include "ita/variation/rov.hpp"	     // for RoV, gen_rov
#include "ita/variation/sne.hpp"	     // for sne

namespace ita::genomics
{
using pool_t = meza::pool::pool<qt::u8, qt::u32>;

void find_inversions_new(const bd::VG &g,
			 const std::set<qt::id_t> &to_call_ref_ids,
			 std::map<qt::u32, std::vector<ia::inv_slice>> &res)
{
	auto foo = [](const liteseq::ref_walk *rw,
		      const std::vector<qt::u32> &positions, qt::u32 u,
		      qt::u32 v) -> qt::op_t<qt::u32>
	{
		if (rw->strands[positions[u]] == liteseq::strand::STRAND_FWD)
			return {positions[u], positions[v]};
		else
			return {positions[v], positions[u]};
	};
	using inv_it =
		dynamo::dynamic_interval_tree<qt::u32, qt::op_t<qt::u32>>;

	std::map<qt::u32, inv_it> hap_idx_to_it;
	for (qt::u32 hap_idx : to_call_ref_ids) {
		std::vector<ia::inv_slice> &inv_slices = res[hap_idx];
		dynamo::dynamic_interval_tree<qt::u32, qt::op_t<qt::u32>> t;
		const liteseq::ref_walk *rw = g.get_ref_vec(hap_idx)->walk;

		for (qt::u32 v_idx{}; v_idx < g.vtx_count(); v_idx++) {
			const std::vector<qt::u32> &positions =
				g.get_vertex_ref_idxs(v_idx, hap_idx);

			qt::u32 N = positions.size();

			if (N < 2)
				continue;

			std::optional<qt::op_t<qt::u32>> d{std::nullopt};
			for (qt::u32 i{}; i < N - 1; i++) {
				if (rw->strands[positions[i]] !=
				    rw->strands[positions[i + 1]]) {
					d = qt::op_t<qt::u32>{i, i + 1};
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
		std::vector<std::pair<qt::u32, qt::u32>> r =
			t.covered_regions();

		for (auto [s, e] : r) {
			qt::u32 len = e - s + 1;
			qt::slice fwd_slice{s, len};
			auto [_, s_rev] = t.values_at(s).front();
			qt::u32 rev_start = s_rev - len + 1;
			qt::slice_t rev_slice{rev_start, len};

			inv_slices.emplace_back(
				ia::inv_slice{rw, hap_idx, fwd_slice.start(),
					      rev_slice.start(), len});
		}
	}
}

void comp_expeditions(const bd::VG &g, std::vector<ir::RoV> &all_rovs,
		      qt::idx_t start, qt::idx_t count,
		      const std::set<qt::id_t> &to_call_ref_ids, pool_t &p,
		      ita::at_matrix::rov_job_batch &batch,
		      std::vector<ia::trek> &treks)
{
	const std::size_t N = all_rovs.size();
	for (qt::idx_t i{start}; i < start + count && i < N; i++) {
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
		     const std::set<qt::id_t> &to_call_ref_ids,
		     bq::bounded_queue<iv::VcfRecIdx> &q,
		     const core::config &app_config)
{
	if (app_config.verbosity() > 0)
		log_info("Generating regions of variation (RoVs)");

	// Parse genomic region if specified
	std::optional<ir::genomic_region> region = std::nullopt;
	if (app_config.has_genomic_region()) {
		region = ir::parse_genomic_region(
			app_config.get_genomic_region().value());
		if (!region.has_value()) {
			log_fatal("Failed to parse genomic region");
			std::exit(EXIT_FAILURE);
		}
	}

	std::vector<ir::RoV> all_rovs =
		ir::gen_rov(pvsts, g, to_call_ref_ids, region);

	const std::size_t CHUNK_SIZE = app_config.get_chunk_size();
	const std::size_t N = all_rovs.size();
	qt::u32 total_chunks = (N + CHUNK_SIZE - 1) / CHUNK_SIZE;

	if (app_config.verbosity() > 0)
		log_info("Processing chunks. Chunk count: {}", total_chunks);

	ise::pin_cushion pc;
	std::vector<ia::trek> treks;
	std::vector<ist::st> i_trees;
	std::map<qt::u32, std::vector<ia::inv_slice>> inv_slices;

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
		for (qt::idx_t base{}; base < N; base += CHUNK_SIZE) {
			qt::u32 end = std::min(base + CHUNK_SIZE, N);
			qt::u32 count = end - base;
			qt::u32 chunk_num = (base / CHUNK_SIZE) + 1;

			if (app_config.verbosity() > 0)
				log_info("\t%ul/%ul", chunk_num, total_chunks);

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
				log_info("Processing inversions");

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
