#include "./call.hpp"


namespace povu::subcommands::call {

/**
 * loop through the .pvst files and read them
 */
void read_pvsts(const core::config &app_config, std::vector<pvst::Tree> &pvsts) {
  std::string fn_name = pv_cmp::format("[povu::subcommands::{}]", __func__);

  // get the list of files in the forest dir that end in .pvst
  std::vector<fs::path> fps = pic::get_files(app_config.get_forest_dir(), ".pvst");

  if (fps.empty()) {
    std::cerr << fn_name << " Could not find pvst files in " << app_config.get_forest_dir() << std::endl;
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

pt::status_t get_refs(core::config &app_config) {
  if (app_config.get_refs_input_fmt() == core::input_format_e::file_path) {
    std::vector<std::string> refs;
    pic::read_lines_to_vec_str(app_config.get_references_txt(), &refs);
    app_config.set_reference_paths(std::move(refs));
    return 0;
  }
  else if (app_config.get_refs_input_fmt() == core::input_format_e::params) {
    // If path prefixes are provided, we need to read all paths from GFA and filter by prefix
    if (!app_config.get_path_prefixes().empty()) {
      std::vector<std::string> filtered_refs = filter_paths_by_prefix(app_config);
      app_config.set_reference_paths(std::move(filtered_refs));
    }
    // If explicit reference paths are provided, they're already set
    return 0;
  }

  return -1;
}


void do_call(core::config &app_config) {
  std::string fn_name = pv_cmp::format("[povu::main::{}]", __func__);

  // ----------------------------------------------------
  // parallel read for the graph, flubbles and references
  // ----------------------------------------------------
  bd::VG *g { nullptr };
  std::vector<pvst::Tree> pvsts;
  pt::status_t _;

  // read graph & refs
  std::thread t1([&] {
    get_refs(app_config);
    g = povu::subcommands::common::get_vg(app_config);
  });

  // read PVST
  std::thread t2([&] { read_pvsts(app_config, pvsts); });

  t1.join();
  t2.join();

#ifdef DEBUG
  if (app_config.verbosity() > 1) {
    g->summary(false);
    INFO("flubble count: {}", pvsts.size());
    INFO("reference count: {}", app_config.get_reference_paths().size());
  }
#endif

  pgv::VcfRecIdx vcf_recs = pg::gen_vcf_rec_map(pvsts, *g, app_config.thread_count());
  piv::write_vcfs(vcf_recs, *g, app_config);

  return;
}

} // namespace povu::subcommands
