#include "./call.hpp"
#include <set>
#include <vector>


namespace povu::subcommands::call {

std::vector<pgt::flubble> get_can_flubbles(const core::config &app_config) {
  std::string fn_name = std::format("[povu::subcommands::{}]", __func__);
  // get the list of files in the forest dir that end in .flb
  std::vector<fs::path> flbs = pic::get_files(app_config.get_forest_dir(), ".flb");

  if (flbs.empty()) {
    std::cerr << fn_name
              << " Could not find flb files in " << app_config.get_forest_dir()
              << std::endl;
    exit(EXIT_FAILURE);
  }

  // TODO: [c] parallelise
  std::vector<pgt::flubble> can_flbs; // flubbles in a given file
  for (std::size_t i{}; i < flbs.size(); i++) {
    //std::cerr << std::format("Reading flubble file: {}\n", flbs[i].string());
    std::vector<pgt::flubble> res = povu::io::bub::read_canonical_fl(flbs[i].string());
    can_flbs.insert(can_flbs.end(), res.begin(), res.end());
  }

  return can_flbs;
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

inline std::set<pt::id_t> get_ref_ids(const bd::VG &g, const core::config &app_config) {
  const std::vector<std::string> &refs = app_config.get_reference_paths();
  std::set<pt::id_t> ref_ids;

  for (const auto &[ref_id, ref_name] : g.get_refs()) {
    if (std::find(refs.begin(), refs.end(), ref_name) != refs.end()) {
      ref_ids.insert(ref_id);
    }
  }
  return ref_ids;
}

void do_call(core::config &app_config) {
  std::string fn_name = std::format("[povu::main::{}]", __func__);

  /* parallel read for the graph, flubbles and references */

  bd::VG *g { nullptr };
  std::vector<pgt::flubble> canonical_flubbles;
  pt::status_t _;

  std::thread t1([&] {
    get_refs(app_config);
    g = pcs::get_vg(app_config);
  });
  std::thread t2([&] { canonical_flubbles = get_can_flubbles(app_config); });
  //std::thread t3([&] { _ = get_refs(app_config); });

  t1.join();
  t2.join();
  //t3.join();

  if (true) { // debug
    g->summary();
    std::cerr << "flubble count = " << canonical_flubbles.size() << "\n";
    std::cerr << "reference count = " << app_config.get_reference_paths().size() << "\n";
  }

  std::set<pt::id_t> ref_ids = get_ref_ids(*g, app_config);

  pvt::VcfRecIdx vcf_recs = pg::gen_vcf_rec_map(canonical_flubbles, *g, app_config);

  piv::write_vcfs(vcf_recs, *g, app_config);

  return;
}

} // namespace povu::subcommands
