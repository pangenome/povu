#include "ita/genomics/genomics.hpp"

#include <algorithm> // for min, max
#include <cmath>     // for ceil
#include <cstddef>   // for size_t
#include <cstdlib>   // for std::max, exit, EXIT_FAILURE
#include <optional>
#include <set>	   // for set
#include <utility> // for move
#include <vector>  // for vector

#include "dynamo/dynamo.hpp"		     // for dynamic_interval_tree
#include "ita/convolutions/convolutions.hpp" // for run_convs
#include "ita/genomics/allele.hpp"	     // for Exp, comp_itineraries
#include "ita/genomics/vcf.hpp"		     // for VcfRecIdx, gen_vcf_records
#include "ita/variation/overlay.hpp"	     // for comp_itineraries3, sub_inv
#include "ita/variation/rov.hpp"	     // for RoV, gen_rov
#include "ita/variation/sne.hpp"	     // for sne
#include "povu/common/app.hpp"		     // for config
#include "povu/common/core.hpp"		     // for pt, idx_t, id_t
#include "povu/common/log.hpp"		     // for ERR
#include "povu/common/thread.hpp"	     // for thread_pool, task_group
#include "povu/refs/refs.hpp"		     // for lq_strand_to_char

namespace ita::genomics
{
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
		      const std::set<pt::id_t> &to_call_ref_ids,
		      std::vector<ia::trek> &treks)
{
	const std::size_t N = all_rovs.size();
	for (pt::idx_t i = start; i < start + count && i < N; ++i) {
		ir::RoV &rov = all_rovs[i];

		std::optional<ia::trek> opt_tk =
			ita::convolutions::comp_expedition(g, rov,
							   to_call_ref_ids);

		if (opt_tk.has_value())
			treks.emplace_back(std::move(opt_tk.value()));
	}

	return;
}

void comp_expeditions_serial(const bd::VG &g, std::vector<ir::RoV> &all_rovs,
			     pt::idx_t start, pt::idx_t count,
			     const std::set<pt::id_t> &to_call_ref_ids,
			     ise::pin_cushion &pc, std::vector<ia::trek> &treks)
{
	const std::size_t N = all_rovs.size();
	for (pt::idx_t i = start; i < start + count && i < N; ++i) {
		ir::RoV &rov = all_rovs[i];
		auto rov_treks =
			po::overlay_generic(g, rov, to_call_ref_ids, pc);

		for (auto &tk : rov_treks)
			treks.emplace_back(std::move(tk));
	}

	return;
}

struct ThreadSplit {
	std::size_t outer;
	std::size_t inner;
};

// Prefer inner >> outer; cap outer to a small number (default 2).
inline ThreadSplit
split_threads(std::size_t total, std::size_t outer_cap = 2,
	      double outer_ratio = 0.25) // use ≤25% of pool for outer
{
	total = std::max<std::size_t>(1, total);
	if (total == 1)
		return {1, 0}; // outer runs; inner is serial

	// Suggested outer from ratio, but no more than outer_cap
	std::size_t suggested_outer = std::max<std::size_t>(
		1, static_cast<std::size_t>(std::ceil(total * outer_ratio)));
	std::size_t outer = std::min(outer_cap, suggested_outer);

	// Always leave at least one thread for inner
	if (outer >= total)
		outer = total - 1;

	std::size_t inner = total - outer; // guaranteed ≥ 1
	return {outer, inner};
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

	// set up thread pool & decide on thread split
	std::size_t thread_count = app_config.thread_count();
	povu::thread::thread_pool pool(thread_count);
	// auto [outer, inner] = split_threads(pool.size());

	ise::pin_cushion pc;
	std::vector<ia::trek> treks;
	std::vector<ist::st> i_trees;
	// std::map<pt::u32, ita::interval_tree::interval_tree> invs;
	std::map<pt::u32, std::vector<ia::inv_slice>> inv_slices;

	try {
		for (pt::idx_t base{}; base < N; base += CHUNK_SIZE) {
			pt::u32 end = std::min(base + CHUNK_SIZE, N);
			pt::u32 count = end - base;
			pt::u32 chunk_num = (base / CHUNK_SIZE) + 1;

			if (app_config.verbosity() > 0)
				INFO("\t{}/{}", chunk_num, total_chunks);

			// comp_expeditions_serial(g, all_rovs, base, count,
			//			to_call_ref_ids, pc, treks);

			comp_expeditions(g, all_rovs, base, count,
					 to_call_ref_ids, treks);

			// if (chunk_num == CHUNK_SIZE)
			//	i_trees = ise::sne(g, pc, to_call_ref_ids);

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
