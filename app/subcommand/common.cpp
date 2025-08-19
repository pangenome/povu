#include "./common.hpp"

namespace  povu::subcommands::common {
/**
 * Read the input gfa into a bidirected variation graph
 */
bd::VG *get_vg(const core::config &app_config) {
  std::chrono::duration<double> timeRefRead;
  auto t0 = pt::Time::now();
  bd::VG *g = povu::io::from_gfa::to_bd(app_config);

#ifdef DEBUG
  if (app_config.verbosity() > 1) {
    timeRefRead = pt::Time::now() - t0;
    INFO("Read GFA file {} into a bidirected graph", app_config.get_input_gfa());
    INFO("Time spent reading the GFA {:.2f} sec", timeRefRead.count());
    t0 = pt::Time::now();
  }
#endif

  return g;
}
}
