#include "./call.hpp"


namespace povu::subcommands::call {

/**
 * loop through the .pvst files and read them
 */
void read_pvsts(const core::config &app_config, std::vector<pvtr::Tree> &pvsts) {
  std::string fn_name = std::format("[povu::subcommands::{}]", __func__);

  // get the list of files in the forest dir that end in .pvst
  std::vector<fs::path> fps = pic::get_files(app_config.get_forest_dir(), ".pvst");

  if (fps.empty()) {
    std::cerr << fn_name << " Could not find pvst files in " << app_config.get_forest_dir() << std::endl;
    exit(EXIT_FAILURE);
  }

  // TODO: [c] parallelise
  // loop through the .pvst files and read them
  for (std::size_t i{}; i < fps.size(); i++) {
    pvtr::Tree pvst = povu::io::pvst::read_pvst(fps[i].string());
    pvst.comp_heights();
    pvsts.push_back(std::move(pvst));
  }
}

pt::status_t get_refs(core::config &app_config) {
  if (app_config.get_refs_input_fmt() != core::input_format_e::file_path) {
    return -1;
  }

  std::vector<std::string> refs;
  pic::read_lines_to_vec_str(app_config.get_references_txt(), &refs);
  app_config.set_reference_paths(std::move(refs));

  return 0;
}


void do_call(core::config &app_config) {
  std::string fn_name = std::format("[povu::main::{}]", __func__);

  // ----------------------------------------------------
  // parallel read for the graph, flubbles and references
  // ----------------------------------------------------
  bd::VG *g { nullptr };
  std::vector<pvtr::Tree> pvsts;
  pt::status_t _;

  // read graph & refs
  std::thread t1([&] {
    get_refs(app_config);
    g = pcs::get_vg(app_config);
  });

  // read PVST
  std::thread t2([&] { read_pvsts(app_config, pvsts); });

  t1.join();
  t2.join();

#ifdef DEBUG
  if (true) {
    g->summary(false);
    std::cerr << "flubble count = " << pvsts.size() << "\n";
    std::cerr << "reference count = " << app_config.get_reference_paths().size() << "\n";
  }
#endif

  //std::set<pt::id_t> ref_ids = get_ref_ids(*g, app_config);
  pvt::VcfRecIdx vcf_recs = pg::gen_vcf_rec_map(pvsts, *g);
  piv::write_vcfs(vcf_recs, *g, app_config);

  return;
}

} // namespace povu::subcommands
