#include "ita/genomics/genomics.hpp"

#include <algorithm> // for min, max
#include <cmath>     // for ceil
#include <cstddef>   // for size_t
#include <cstdlib>   // for std::max, exit, EXIT_FAILURE
#include <numeric>   // for accumulate
#include <set>	     // for set
#include <utility>   // for move
#include <vector>    // for vector

#include "ita/genomics/allele.hpp"     // for Exp, comp_itineraries
#include "ita/genomics/vcf.hpp"	       // for VcfRecIdx, gen_vcf_records
#include "ita/graph/interval_tree.hpp" // for interval_tree
#include "ita/variation/overlay.hpp"   // for comp_itineraries3, sub_inv
#include "ita/variation/rov.hpp"       // for RoV, gen_rov
#include "ita/variation/sne.hpp"       // for sne
#include "povu/common/app.hpp"	       // for config
#include "povu/common/core.hpp"	       // for pt, idx_t, id_t
#include "povu/common/log.hpp"	       // for ERR
#include "povu/common/thread.hpp"      // for thread_pool, task_group

namespace ita::genomics
{
void find_inversions(
	const bd::VG &g, const std::set<pt::id_t> &to_call_ref_ids,
	std::map<pt::u32, ita::interval_tree::interval_tree> &hap_idx_to_it)
{
	for (pt::u32 hap_idx : to_call_ref_ids) {
		// pr::Ref r = g.get_ref_by_id(hap_idx);
		// std::cerr << r.tag() << "\n";
		// std::cerr << "[";
		ita::interval_tree::interval_tree it{hap_idx};
		for (pt::u32 v_idx{}; v_idx < g.vtx_count(); v_idx++) {
			const std::vector<pt::u32> &positions =
				g.get_vertex_ref_idxs(v_idx, hap_idx);

			pt::u32 N = positions.size();

			if (N < 2)
				continue;

			const liteseq::ref_walk *rw =
				g.get_ref_vec(hap_idx)->walk;

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

			auto [f, s] = (rw->strands[positions[u]] ==
				       liteseq::strand::STRAND_FWD)
					      ? pt::op_t<pt::u32>{positions[u],
								  positions[v]}
					      : pt::op_t<pt::u32>{positions[v],
								  positions[u]};

			it.add(f, s, 1);

			// pt::u32 f_vid = rw->v_ids[f];
			// pt::u32 s_vid = rw->v_ids[s];

			// if (f_vid != s_vid)
			//	continue;

			// liteseq::strand f_strand = rw->strands[f];
			// liteseq::strand s_strand = rw->strands[s];

			// if (f_strand == s_strand)
			//	continue;

			// std::cerr << "{" << g.v_idx_to_id(v_idx) << " " <<
			// s_vid
			//	  << " " << f_vid << "}, ";

			// auto [fwd_idx, rev_idx] =
			//	(f_strand == liteseq::strand::STRAND_FWD)
			//		? pt::op_t<pt::u32>{f, s}
			//		: pt::op_t<pt::u32>{s, f};

			// invs.add_inv_slice(v_idx, hap_idx,
			//		   ia::inv_slice{rw, hap_idx, fwd_idx,
			//				 rev_idx, 1});
		}

		hap_idx_to_it.emplace(hap_idx, std::move(it));
		// std::cerr << "]\n";
	}

	// for (pt::u32 hap_idx{}; hap_idx < g.get_hap_count(); hap_idx++) {
	//	const liteseq::ref_walk *rw = g.get_ref_vec(hap_idx)->walk;
	//	pr::Ref r = g.get_ref_by_id(hap_idx);
	//	std::cerr << r.tag() << "\n";
	//	std::cerr << "[";
	//	for (const auto &[v_idx, inv_slices] :
	//	     invs.get_slices_by_hap_idx(hap_idx)) {
	//		std::cerr << g.v_idx_to_id(v_idx) << ", ";
	//	}
	//	std::cerr << "]\n";
	// }

	// return hap_idx_to_it;
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

	if (app_config.verbosity() > 0)
		INFO("No. of chunks: {}", N);

	// set up thread pool & decide on thread split
	std::size_t thread_count = app_config.thread_count();
	povu::thread::thread_pool pool(thread_count);
	// auto [outer, inner] = split_threads(pool.size());

	ise::pin_cushion pc;
	std::vector<ia::trek> treks;

	std::vector<ist::st> i_trees;

	std::map<pt::u32, ita::interval_tree::interval_tree> invs;

	try {
		for (pt::idx_t base{}; base < N; base += CHUNK_SIZE) {
			pt::u32 end = std::min(base + CHUNK_SIZE, N);
			pt::u32 count = end - base;
			pt::u32 chunk_num = (base / CHUNK_SIZE) + 1;

			if (app_config.verbosity() > 0)
				INFO("Processing RoV Chunk ({}/{})", chunk_num,
				     (N + CHUNK_SIZE - 1) / CHUNK_SIZE);

			comp_expeditions_serial(g, all_rovs, base, count,
						to_call_ref_ids, pc, treks);

			if (chunk_num == CHUNK_SIZE)
				i_trees = ise::sne(g, pc, to_call_ref_ids);

			iv::VcfRecIdx rs =
				iv::gen_vcf_records(g, treks, i_trees, invs);

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

			find_inversions(g, to_call_ref_ids, invs);

			iv::VcfRecIdx rs =
				iv::gen_vcf_records(g, treks, i_trees, invs);
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
