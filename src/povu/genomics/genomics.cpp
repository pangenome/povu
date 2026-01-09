#include "povu/genomics/genomics.hpp"

#include <algorithm> // for min, max
#include <cmath>     // for ceil
#include <cstddef>   // for size_t
#include <cstdlib>   // for std::max, exit, EXIT_FAILURE
#include <utility>   // for move
#include <vector>

#include "povu/common/app.hpp" // for config
#include "povu/common/core.hpp"
#include "povu/common/log.hpp"	    // for ERR
#include "povu/common/thread.hpp"   // for thread_pool, task_group
#include "povu/genomics/allele.hpp" // for Exp, comp_itineraries
#include "povu/genomics/vcf.hpp"    // for VcfRecIdx, gen_vcf_records
#include "povu/overlay/overlay.hpp" // for comp_itineraries3, sub_inv
#include "povu/overlay/sne.hpp"	    // for sne
#include "povu/variation/rov.hpp"   // for RoV, gen_rov

namespace povu::genomics
{
namespace pga = povu::genomics::allele;
namespace pgv = povu::genomics::vcf;
namespace pvst = povu::pvst;

void comp_expeditions_serial(const bd::VG &g, std::vector<pvr::RoV> &all_rovs,
			     pt::idx_t start, pt::idx_t count,
			     const std::set<pt::id_t> &to_call_ref_ids,
			     pos::pin_cushion &pc,
			     std::vector<pga::trek> &treks)
{
	const std::size_t N = all_rovs.size();
	for (pt::idx_t i = start; i < start + count && i < N; ++i) {
		pvr::RoV &rov = all_rovs[i];
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
		     pbq::bounded_queue<pgv::VcfRecIdx> &q,
		     const core::config &app_config)
{
	// bool prog = app_config.show_progress();

	// Parse genomic region if specified
	std::optional<pvr::genomic_region> region = std::nullopt;
	if (app_config.has_genomic_region()) {
		region = pvr::parse_genomic_region(
			app_config.get_genomic_region().value());
		if (!region.has_value()) {
			ERR("Failed to parse genomic region");
			std::exit(EXIT_FAILURE);
		}
	}

	std::vector<pvr::RoV> all_rovs =
		pvr::gen_rov(pvsts, g, to_call_ref_ids, region);

	// std::cerr << "found " << all_rovs.size();

	const std::size_t CHUNK_SIZE = app_config.get_chunk_size();
	const std::size_t N = all_rovs.size();
	// const std::size_t CHUNK_COUNT = (N + CHUNK_SIZE - 1) / CHUNK_SIZE;

	// set up progress bars
	// ProgressBar chunks_prog_bar{option::Stream{std::cerr}};
	// set_progress_bar_common_opts(&chunks_prog_bar, CHUNK_COUNT);
	// std::string prog_msg; // setup buffer for progress bar messages
	// prog_msg.reserve(128);

	// set up thread pool & decide on thread split
	std::size_t thread_count = app_config.thread_count();
	povu::thread::thread_pool pool(thread_count);
	// auto [outer, inner] = split_threads(pool.size());

	pos::pin_cushion pc;
	std::vector<pga::trek> treks;

	// std::vector<pga::Exp> exps;
	// treks.reserve(CHUNK_SIZE);

	try {
		for (pt::idx_t base{}; base < N; base += CHUNK_SIZE) {
			pt::u32 end = std::min(base + CHUNK_SIZE, N);
			pt::u32 count = end - base;
			pt::u32 chunk_num = (base / CHUNK_SIZE) + 1;

			// if (prog) {
			//	pt::idx_t chunk_num = (base /
			// CHUNK_SIZE) + 1;

			//	prog_msg = pv_cmp::format(
			//		"Processing RoV Chunk ({}/{})",
			//		chunk_num, CHUNK_COUNT);
			//	chunks_prog_bar.set_option(
			//		option::PostfixText{prog_msg});
			//	chunks_prog_bar.set_progress(
			//		static_cast<size_t>(chunk_num));
			// }

			// std::cerr << "A\n";

			comp_expeditions_serial(g, all_rovs, base, count,
						to_call_ref_ids, pc, treks);

			std::vector<poi::it> i_trees;
			if (chunk_num == CHUNK_SIZE)
				i_trees = pos::sne(g, pc, to_call_ref_ids);

			if (treks.empty() && i_trees.empty())
				continue;

			pgv::VcfRecIdx rs =
				pgv::gen_vcf_records(g, treks, i_trees);

			if (!q.push(std::move(rs)))
				break; // queue was closed early

			treks.clear();
		}
	}
	catch (...) {
		q.close(); // make sure consumers wake up on errors
		throw;
	}

	q.close(); // we're done
}
} // namespace povu::genomics
