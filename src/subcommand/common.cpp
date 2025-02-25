#include "./common.hpp"

namespace  povu::subcommands::common {
/**
 * Read the input gfa into a bidirected variation graph
 */
bd::VG *get_vg(const core::config &app_config) {
  std::string fn_name = std::format("[povu::deconstruct::{}]", __func__);
  std::size_t ll = app_config.verbosity(); // log level

  std::chrono::duration<double> timeRefRead;
  auto t0 = pt::Time::now();

  // if (ll > 2) {
  //   std::cerr << std::format("{} Reading graph\n", fn_name);
  // }
  bd::VG *g = povu::io::from_gfa::to_bd(app_config);

  if (ll > 1) {
    timeRefRead = pt::Time::now() - t0;
    povu::utils::report_time(std::cerr, fn_name, "read_gfa", timeRefRead);
    t0 = pt::Time::now();
  }

  return g;
}
}
