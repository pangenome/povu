#include "./call.hpp"
#include "indicators/indeterminate_progress_bar.hpp"
#include "indicators/progress_bar.hpp"
#include <string>
#include <thread>



namespace povu::subcommands::call {

/**
 * loop through the .pvst files and read them
 */
void read_pvsts(const core::config &app_config, std::vector<pvst::Tree> &pvsts) {
  //std::string fn_name = pv_cmp::format("[povu::subcommands::{}]", __func__);

  // get the list of files in the forest dir that end in .pvst
  std::vector<fs::path> fps = pic::get_files(app_config.get_forest_dir(), ".pvst");

  if (fps.empty()) {
    ERR("Could not find pvst files in {}", app_config.get_forest_dir().string());
    exit(EXIT_FAILURE);
  }

  // TODO: [c] parallelise
  // loop through the .pvst files and read them
  for (std::size_t i{}; i < fps.size(); i++) {
    pvst::Tree pvst = povu::io::from_pvst::read_pvst(fps[i].string());
    pvst.comp_heights();
    pvsts.push_back(std::move(pvst));
  }
}

std::vector<std::string> filter_paths_by_prefix(const core::config &app_config) {
  std::string fn_name = pv_cmp::format("[povu::subcommands::{}]", __func__);

  // Read GFA file and extract all P-line names
  std::ifstream gfa_file(app_config.get_input_gfa());
  if (!gfa_file.is_open()) {
    std::cerr << fn_name << " Could not open GFA file: " << app_config.get_input_gfa() << std::endl;
    exit(EXIT_FAILURE);
  }

  std::vector<std::string> all_paths;
  std::string line;

  while (std::getline(gfa_file, line)) {
    if (line.empty() || line[0] != 'P') continue;

    // P-line format: P<tab>path_name<tab>path_in_graph<tab>cigars
    std::istringstream ss(line);
    std::string token;

    // Skip the 'P' token
    std::getline(ss, token, '\t');
    if (token != "P") continue;

    // Get path name
    if (std::getline(ss, token, '\t')) {
      all_paths.push_back(token);
    }
  }
  gfa_file.close();

  // Filter paths by prefixes
  std::vector<std::string> filtered_paths;
  const auto& prefixes = app_config.get_path_prefixes();

  for (const std::string& path : all_paths) {
    for (const std::string& prefix : prefixes) {
      if (path.substr(0, prefix.length()) == prefix) {
        filtered_paths.push_back(path);
        break; // Found matching prefix, no need to check others
      }
    }
  }

  if (filtered_paths.empty()) {
    std::cerr << fn_name << " Warning: No paths found matching the provided prefixes" << std::endl;
  }

  return filtered_paths;
}

void get_refs(core::config &app_config) {
  if (app_config.get_refs_input_fmt() == core::input_format_e::file_path) {
    std::vector<std::string> refs;
    pic::read_lines_to_vec_str(app_config.get_references_txt(), &refs);
    app_config.set_reference_paths(std::move(refs));
    return;
  }
  else if (app_config.get_refs_input_fmt() == core::input_format_e::params) {
    // If path prefixes are provided, we need to read all paths from GFA and filter by prefix
    if (!app_config.get_path_prefixes().empty()) {
      std::vector<std::string> filtered_refs = filter_paths_by_prefix(app_config);
      app_config.set_reference_paths(std::move(filtered_refs));
    }
    // If explicit reference paths are provided, they're already set
    return;
  }


  ERR("Could not set refs");
  std::exit(EXIT_FAILURE);
}




void do_call(core::config &app_config) {

  // ----------------------------------------------------
  // parallel read for the graph, flubbles and references
  // ----------------------------------------------------
  bd::VG *g { nullptr };
  std::vector<pvst::Tree> pvsts;

  // read graph
  std::thread read_graph([&] {
    g = povu::io::from_gfa::to_bd(app_config);
  });

  // read PVST
  std::thread read_pvsts_async([&] {
    read_pvsts(app_config, pvsts);
  });

  // read references either from file or directly from params
  std::thread get_refs_async([&] {
    get_refs(app_config);
  });

  // create VCF output object
  // assumes no duplicates in samples
  read_graph.join();
  get_refs_async.join();
  const std::vector<std::string> &samples = app_config.get_reference_paths();

#ifdef DEBUG
  assert(samples.size() > 0 && "No reference samples found");
#endif

  std::set<pt::id_t> vcf_ref_ids; // ref IDs that we need to output VCF for
  for (const auto &sample : samples) {
    std::set<pt::id_t> ref_ids = g->get_refs_in_sample(sample);
    vcf_ref_ids.insert(ref_ids.begin(), ref_ids.end());
  }

  piv::VcfOutput vout = app_config.get_stdout_vcf()
    ? piv::VcfOutput::to_stdout()
    : piv::VcfOutput::to_split_files(app_config.get_reference_paths(), std::string(app_config.get_output_dir()));

  std::thread init_vcfs_async([&] {
    piv::init_vcfs(*g, samples, vout);
  });

  read_pvsts_async.join();

  DynamicProgress<ProgressBar> prog_bars;
  ProgressBar bar1;
  set_progress_bar_common_opts(&bar1);
  std::size_t bar_idx = prog_bars.push_back(bar1);

  // if running out of memory, reduce the capacity and/or the chunk size
  const std::size_t QUEUE_CAPACITY = app_config.get_queue_len();
  pbq::bounded_queue<pgv::VcfRecIdx> q(QUEUE_CAPACITY);

  // start producer in its own thread
  std::thread producer([&] {
    try {
      pg::gen_vcf_rec_map(pvsts, *g, q, prog_bars, bar_idx, app_config.thread_count(), app_config);
    } catch (...) {
      q.close(); // make sure consumers wake up on errors
      throw;
    }
  });


  // // add dynamic progress bar for VCF writing
  // DynamicProgress<IndeterminateProgressBar> ind_bars;
  // IndeterminateProgressBar bar2;
  // set_progress_bar_ind(&bar2);
  // bar2.set_option(indicators::option::PostfixText{"Genereating VCFs..."});
  // std::size_t ind_bar_idx = ind_bars.push_back(bar2);
  // std::thread prog([&] {
  //   while (!q.closed()) {
  //     ind_bars[ind_bar_idx].tick();
  //     std::this_thread::sleep_for(std::chrono::milliseconds(100));
  //   }
  //   if (q.closed()) {
  //     bar2.mark_as_completed();
  //   }
  // });
  // prog.detach();

  // consumer on this thread
  while (auto opt_rec_idx = q.pop()) {
    piv::write_vcfs(*opt_rec_idx, *g, vcf_ref_ids, vout, app_config);
    // sleep for a bit to give more threads to producer
    //std::this_thread::sleep_for(std::chrono::milliseconds(10));
  }

  // make sure VCF are initialised before producer finishes
  init_vcfs_async.join();

  // wait for producer to finish
  producer.join();

  // just in case
  vout.flush_all();

  return;
}

} // namespace povu::subcommands
