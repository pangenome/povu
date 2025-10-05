#include "./call.hpp"

#include <cassert>    // for assert
#include <cstdlib>    // for size_t, exit, EXIT_FAILURE
#include <filesystem> // for path
#include <optional>   // for optional
#include <set>	      // for set
#include <string>     // for basic_string, string, cha...
#include <thread>     // for thread
#include <utility>    // for move
#include <vector>     // for vector

#include "indicators/dynamic_progress.hpp" // for DynamicProgress
#include "indicators/progress_bar.hpp"	   // for ProgressBar
#include "povu/common/bounded_queue.hpp"   // for pbq, bounded_queue
#include "povu/common/core.hpp"		   // for pt, id_t
#include "povu/common/log.hpp"		   // for ERR
#include "povu/common/progress.hpp"	   // for set_progress_bar_common_opts
#include "povu/genomics/genomics.hpp"	   // for gen_vcf_rec_map
#include "povu/genomics/vcf.hpp"	   // for VcfRecIdx
#include "povu/graph/bidirected.hpp"	   // for VG, bd
#include "povu/graph/pvst.hpp"		   // for Tree
#include "povu/io/common.hpp"		   // for get_files, read_lines_to_...
#include "povu/io/from_gfa.hpp"		   // for to_bd
#include "povu/io/from_pvst.hpp"	   // for read_pvst
#include "povu/io/to_vcf.hpp"		   // for VcfOutput, init_vcfs, wri...

namespace povu::subcommands::call
{
using namespace povu::progress;
namespace fs = std::filesystem;
namespace pgv = povu::genomics::vcf;
namespace pg = povu::genomics;
namespace pic = povu::io::common;
namespace piv = povu::io::to_vcf;

/**
 * loop through the .pvst files and read them
 */
void read_pvsts(const core::config &app_config, std::vector<pvst::Tree> &pvsts)
{

	// get the list of files in the forest dir that end in .pvst
	std::vector<fs::path> fps =
		pic::get_files(app_config.get_forest_dir(), ".pvst");

	if (fps.empty()) {
		ERR("Could not find pvst files in {}",
		    app_config.get_forest_dir().string());
		exit(EXIT_FAILURE);
	}

	// TODO: [c] parallelise
	// loop through the .pvst files and read them
	for (std::size_t i{}; i < fps.size(); i++) {
		pvst::Tree pvst =
			povu::io::from_pvst::read_pvst(fps[i].string());
		pvst.comp_heights();
		pvsts.push_back(std::move(pvst));
	}
}

void get_ref_prefixes_from_file(core::config &app_config)
{
	if (app_config.get_refs_input_fmt() ==
	    core::input_format_e::file_path) {
		std::vector<std::string> refs;
		pic::read_lines_to_vec_str(app_config.get_references_txt(),
					   &refs);
		for (auto &&r : refs) {
			app_config.add_ref_name_prefix(r);
		}
	}
}

void do_call(core::config &app_config)
{

	// ----------------------------------------------------
	// parallel read for the graph, flubbles and references
	// ----------------------------------------------------
	bd::VG *g{nullptr};
	std::vector<pvst::Tree> pvsts;

	// read graph
	std::thread read_graph([&]
			       { g = povu::io::from_gfa::to_bd(app_config); });

	// read PVST
	std::thread read_pvsts_async([&] { read_pvsts(app_config, pvsts); });

	// read references either from file or directly from params
	std::thread get_refs_async([&]
				   { get_ref_prefixes_from_file(app_config); });

	// create VCF output object
	// assumes no duplicates in samples
	read_graph.join();
	get_refs_async.join();
	const std::vector<std::string> &samples =
		app_config.get_ref_name_prefixes();

#ifdef DEBUG
	assert(samples.size() > 0 && "No reference samples found");
#endif

	std::set<pt::id_t>
		vcf_ref_ids; // ref IDs that we need to output VCF for
	for (std::string_view sample : samples) {
		std::set<pt::id_t> ref_ids = g->get_refs_in_sample(sample);
		vcf_ref_ids.insert(ref_ids.begin(), ref_ids.end());
	}

#ifdef DEBUG
	assert(vcf_ref_ids.size() > 0 &&
	       "No reference IDs found for VCF output");
#endif

	piv::VcfOutput vout =
		app_config.get_stdout_vcf()
			? piv::VcfOutput::to_stdout()
			: piv::VcfOutput::to_split_files(
				  app_config.get_ref_name_prefixes(),
				  std::string(app_config.get_output_dir()));

	std::thread init_vcfs_async([&] { piv::init_vcfs(*g, samples, vout); });

	read_pvsts_async.join(); // make sure pvsts are read

	DynamicProgress<ProgressBar> bars;
	ProgressBar chunks_prog_bar;
	set_progress_bar_common_opts(&chunks_prog_bar);
	std::size_t chunks_prog_bar_idx = bars.push_back(chunks_prog_bar);

	// if running out of memory, reduce the capacity and/or the chunk size
	const std::size_t QUEUE_CAPACITY = app_config.get_queue_len();
	pbq::bounded_queue<pgv::VcfRecIdx> q(QUEUE_CAPACITY);

	// start producer in its own thread
	std::thread producer(
		[&]
		{
			try {
				pg::gen_vcf_rec_map(pvsts, *g, vcf_ref_ids, q,
						    bars, chunks_prog_bar_idx,
						    app_config);
			}
			catch (...) {
				q.close(); // make sure consumers wake up on
					   // errors
				throw;
			}
		});

	// consumer on this thread
	while (auto opt_rec_idx = q.pop()) {
		piv::write_vcfs(*opt_rec_idx, *g, vcf_ref_ids, vout,
				app_config);
	}

	// make sure VCF are initialised before producer finishes
	init_vcfs_async.join();

	// wait for producer to finish
	producer.join();

	// just in case
	vout.flush_all();

	delete g;
	g = nullptr;

	return;
}

} // namespace povu::subcommands::call
