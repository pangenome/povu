#include "./decompose.hpp"

#include <cstddef>    // for size_t
#include <cstdint>    // for uint32_t
#include <functional> // for ref, reference_wrapper
#include <iostream>   // for basic_ostream, cerr, operat...
#include <string>     // for basic_string, operator<<
#include <thread>     // for thread
#include <utility>    // for get, make_pair, pair
#include <vector>     // for vector

#include "fmt/core.h"			 // for format
#include "povu/algorithms/concealed.hpp" // for find_concealed
#include "povu/algorithms/flubbles.hpp"	 // for flubbles
#include "povu/algorithms/midi.hpp"	 // for find_midi
#include "povu/algorithms/parallel.hpp"	 // for find_parallel
#include "povu/algorithms/smothered.hpp" // for find_smothered
#include "povu/algorithms/tiny.hpp"	 // for find_tiny
#include "povu/common/app.hpp"		 // for config
#include "povu/common/compat.hpp"	 // for pv_cmp, format
#include "povu/graph/bidirected.hpp"	 // for bidirected
#include "povu/graph/pvst.hpp"		 // for pvst
#include "povu/graph/spanning_tree.hpp"	 // for spanning_tree
#include "povu/graph/tree_utils.hpp"	 // for tree_utils
#include "povu/io/from_gfa.hpp"		 // for to_bd
#include "povu/io/to_gfa.hpp"
#include "povu/io/to_pvst.hpp" // for write_pvst, pv_to_pvst

namespace povu::subcommands::decompose
{
namespace ptu = povu::tree_utils;
namespace pfl = povu::flubbles;

void decompose_component(bd::VG *g, std::size_t component_id,
			 const core::config &app_config)
{
	const std::string fn_name =
		pv_cmp::format("[{}::{}]", MODULE, __func__);

#ifdef DEBUG
	if (app_config.verbosity() > 4) {
		std::cerr << "\n";
		g->print_dot(std::cerr);
		std::cerr << "\n";
	}
#endif

	pst::Tree st = pst::Tree::from_bd(*g);
	delete g;

	ptu::tree_meta tm = ptu::gen_tree_meta(st);

	pvst::Tree flubble_tree = pfl::find_flubbles(st, app_config);

#ifdef DEBUG
	if (app_config.verbosity() > 4) {
		std::cerr << "\n";
		st.print_dot(std::cerr);
		std::cerr << "\n";
	}
#endif

	if (app_config.find_subflubbles()) {
		povu::tiny::find_tiny(st, flubble_tree, tm);
		povu::parallel::find_parallel(st, flubble_tree, tm);
		povu::concealed::find_concealed(st, flubble_tree, tm);
		povu::midi::find_midi(st, flubble_tree, tm);
		povu::smothered::find_smothered(st, flubble_tree, tm);
	}

	pv_to_pvst::write_pvst(flubble_tree, std::to_string(component_id),
			       app_config);

	return;
}

std::pair<uint32_t, uint32_t> thread_count(const core::config &app_config,
					   std::size_t item_count)
{
	/* Divide the number of components into chunks for each thread */
	unsigned int total_threads = std::thread::hardware_concurrency();

	pt::u32 conf_num_threads =
		static_cast<std::size_t>(app_config.thread_count());
	uint32_t num_threads = (conf_num_threads > total_threads)
				       ? total_threads
				       : conf_num_threads;
	uint32_t chunk_size = item_count / num_threads;

	return std::make_pair(num_threads, chunk_size);
}

void do_decompose(const core::config &app_config)
{
	std::size_t ll = app_config.verbosity();	   // ll for log level
	bd::VG *g = povu::io::from_gfa::to_bd(app_config); // read graph

	if (ll > 1)
		INFO("Finding components");

	std::vector<bd::VG *> components = bd::VG::componetize(*g);

	delete g;

	if (ll > 1)
		INFO("Found {} components", components.size());

	std::pair<pt::u32, pt::u32> thread_config =
		thread_count(app_config, components.size());

	pt::u32 num_threads = thread_config.first;
	pt::u32 chunk_size = thread_config.second;

	/* Create and launch threads */
	std::vector<std::thread> threads(num_threads);
	pt::u32 start, end;
	for (pt::u32 thread_idx{}; thread_idx < num_threads; ++thread_idx) {
		start = thread_idx * chunk_size;
		end = (thread_idx == num_threads - 1)
			      ? components.size()
			      : (thread_idx + 1) * chunk_size;

		threads[thread_idx] = std::thread(
			[start, end, num_threads, app_config, ll, &components]
			{
				for (pt::u32 i{start}; i < end; i++) {

					pt::u32 component_id{i + 1};

					if (app_config.verbosity())
						INFO("Handling component: {}",
						     component_id);

					pt::u32 N = components[i]->vtx_count();
					if (N < 3) {
						// clang-format off
						if (ll > 2)
							INFO("Skipping component {} because it is too small. (size: {})", component_id, N);
						// clang-format on
						continue;
					}

					if (ll > 3 && num_threads == 1)
						components[i]->summary(false);

					// pass app_config by reference
					decompose_component(
						components[i], component_id,
						std::ref(app_config));
				}
			});
	}

	// Wait for all threads to finish
	for (auto &thread : threads)
		thread.join();

	return;
}

} // namespace povu::subcommands::decompose
