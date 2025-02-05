#include "./subcommands.hpp"

namespace povu::subcommands {

void do_call(const core::config &app_config) {
  std::string fn_name = std::format("[povu::main::{}]", __func__);

  std::chrono::duration<double> timeRefRead;
  auto t0 = pt::Time::now();

  // -----
  // read the input gfa into a bidirected variation graph
  // -----
  if (app_config.verbosity() > 2)  { std::cerr << std::format ("{} Reading graph\n", fn_name); }

  t0 = pt::Time::now();
  bd::VG bd_vg = io::from_gfa::to_bd(app_config.get_input_gfa().c_str(), app_config);
  if (app_config.verbosity() > 1) {
    timeRefRead = pt::Time::now() - t0;
    povu::utils::report_time(std::cerr, fn_name, "read_gfa2", timeRefRead);
  }

  std::vector<std::filesystem::path> flubble_files = povu::io::generic::get_files(app_config.get_forest_dir(), ".flb");

  std::vector<pgt::flubble> canonical_flubbles;

  // TODO: do in parallel
  for (const auto& fp: flubble_files) {
    std::cerr << std::format("Reading flubble file: {}\n", fp.string());
    auto res = povu::io::bub::read_canonical_fl(fp.string());
    // append the results of the vector with these results
    canonical_flubbles.insert(canonical_flubbles.end(), res.begin(), res.end());
  }

  // ------
  // read from a flubble tree in flb in format
  // -----
  //povu::genomics::call_variants(canonical_flubbles, bd_vg, app_config);

  return;
}

} // namespace povu::subcommands
