#include "povu/genomics/genomics.hpp"

#include <algorithm> // for min, max
// #include <atomic>    // for atomic, memory_order
#include <cmath>   // for ceil
#include <cstddef> // for size_t
#include <cstdlib> // for std::max
// #include <optional>  // for optional, operator==
// #include <string>    // for basic_string, string
#include <utility> // for move
#include <vector>

// #include "fmt/core.h" // for format_to
//  #include "indicators/progress_bar.hpp"
//  #include "indicators/setting.hpp" // for PostfixText
#include "povu/common/app.hpp" // for config
#include "povu/common/core.hpp"
#include "povu/common/thread.hpp"     // for thread_pool, task_group
#include "povu/genomics/allele.hpp"   // for Exp, comp_itineraries
#include "povu/genomics/untangle.hpp" // for untangle_ref_walks
#include "povu/genomics/vcf.hpp"      // for VcfRecIdx, gen_vcf_records
#include "povu/overlay/overlay.hpp"   // for comp_itineraries3, sub_inv
#include "povu/overlay/sne.hpp"	      // for sne
#include "povu/variation/rov.hpp"     // for RoV, gen_rov

namespace povu::genomics
{
// using namespace povu::progress;
namespace pga = povu::genomics::allele;
namespace pgv = povu::genomics::vcf;
namespace put = povu::genomics::untangle;
namespace pvst = povu::pvst;

// /**
//  * Associate walks in an RoV with references
//  */
// pga::Exp exp_frm_rov(const bd::VG &g, const pvr::RoV &rov)
// {
//	pga::Exp exp(&rov);
//	pga::comp_itineraries(g, exp);
//	if (exp.is_tangled())
//		put::untangle_ref_walks(exp);

//	return exp;
// }

std::vector<pga::Exp>
comp_expeditions_serial(const bd::VG &g, const std::vector<pvr::RoV> &all_rovs,
			pt::idx_t start, pt::idx_t count,
			const std::set<pt::id_t> &to_call_ref_ids)
{
	const std::size_t N = all_rovs.size();
	std::vector<pga::Exp> all_exp;
	std::vector<pos::pin_cushion> pcushions;
	all_exp.reserve(N * 2);

	for (pt::idx_t i = start; i < start + count && i < N; ++i) {
		const pvr::RoV &rov = all_rovs[i];
		auto [rov_exps, pcs] =
			po::comp_itineraries3(g, rov, to_call_ref_ids);

		for (auto &si : pcs)
			pcushions.emplace_back(std::move(si));

		for (auto &e : rov_exps) {
			if (e.is_tangled())
				put::untangle_ref_walks(e);

			all_exp.emplace_back(std::move(e));
		}
	}

	if (!pcushions.empty())
		pos::sne(g, pcushions, to_call_ref_ids, all_exp);

	all_exp.shrink_to_fit();

	return all_exp;
}

// std::vector<pga::Exp> comp_expeditions_work_steal(
//	const bd::VG &g, const std::vector<pvr::RoV> &all_rovs, pt::idx_t start,
//	pt::idx_t count, povu::thread::thread_pool &pool,
//	std::size_t outer_concurrency, std::size_t reserve_for_inner)
// {
//	auto base = static_cast<std::size_t>(start);
//	auto want = static_cast<std::size_t>(count);

//	const std::size_t N = count;
//	std::vector<pga::Exp> all_exp(N);

//	// choose K workers; leave some for inner tasks
//	const std::size_t pool_sz = pool.size();
//	if (reserve_for_inner >= pool_sz) {
//		reserve_for_inner = 0;
//	}
//	const std::size_t K = std::max<std::size_t>(
//		1,
//		std::min({outer_concurrency, pool_sz - reserve_for_inner, N}));

//	std::atomic<std::size_t> next{0}; // counts *offsets* within [0, N)
//	const std::size_t block = 64;

//	povu::thread::task_group tg(pool);

//	for (std::size_t k = 0; k < K; ++k) {
//		tg.run(
//			[&, base]
//			{
//				for (;;) {
//					// offset within the chunk
//					const std::size_t off = next.fetch_add(
//						block,
//						std::memory_order_relaxed);
//					if (off >= N) {
//						break;
//					}

//					// global index range [g_begin,
//					// g_end)
//					const std::size_t g_begin = base + off;
//					const std::size_t g_end = std::min(
//						g_begin + block, base + N);

//					// fill results; local index =
//					// global - base
//					for (std::size_t gi = g_begin;
//					     gi < g_end; ++gi) {
//						// const std::size_t li =
//						//	gi - base;
//						// IMPORTANT:
//						// - Pass the RoV
//						// element directly: do
//						// NOT make a local copy
//						// like
//						//   `pvr::RoV rov =
//						//   all_rovs[gi];` and
//						//   then take &rov —
//						//   that would dangle.
//						// - `exp_frm_rov(...)`
//						// returns an Exp
//						// prvalue. Inside that
//						// function,
//						//   the return is
//						//   elided (RVO), and
//						//   here the assignment
//						//   uses
//						//   move-assignment. No
//						//   extra copy is
//						//   performed; no
//						//   `std::move` is
//						//   needed on the RHS.
//						// - Ensure `all_rovs`
//						// outlives all Exp
//						// objects that store
//						// pointers
//						//   into it.

//						// all_exp[li] = exp_frm_rov(
//						//	g, all_rovs[gi]);
//					}
//				}
//			});
//	}

//	tg.wait(); // rethrows first exception if any

//	return all_exp;
// }

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

	std::vector<pvr::RoV> all_rovs =
		pvr::gen_rov(pvsts, g, to_call_ref_ids);

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

	std::vector<pga::Exp> exps;
	exps.reserve(CHUNK_SIZE);

	try {
		for (pt::idx_t base{}; base < N; base += CHUNK_SIZE) {
			const pt::idx_t end = std::min(base + CHUNK_SIZE, N);
			const pt::idx_t count = end - base;

			// if (prog) {
			//	pt::idx_t chunk_num = (base / CHUNK_SIZE) + 1;

			//	prog_msg = pv_cmp::format(
			//		"Processing RoV Chunk ({}/{})",
			//		chunk_num, CHUNK_COUNT);
			//	chunks_prog_bar.set_option(
			//		option::PostfixText{prog_msg});
			//	chunks_prog_bar.set_progress(
			//		static_cast<size_t>(chunk_num));
			// }

			exps = comp_expeditions_serial(g, all_rovs, base, count,
						       to_call_ref_ids);

			if (exps.empty())
				continue;

			pgv::VcfRecIdx rs =
				pgv::gen_vcf_records(g, exps, to_call_ref_ids);

			if (!q.push(std::move(rs)))
				break; // queue was closed early

			exps.clear();
		}
	}
	catch (...) {
		q.close(); // make sure consumers wake up on errors
		throw;
	}

	q.close(); // we're done
}
} // namespace povu::genomics
