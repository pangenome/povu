#include "./subcommands.hpp"

namespace povu::subcommands {

void do_info(const core::config &app_config) {
  std::string fn_name = std::format("[povu::main::{}]", __func__);

  // -----
  // read the input gfa into a bidirected variation graph
  // -----
  bd::VG g = io::from_gfa::to_bd(app_config.get_input_gfa().c_str(), app_config);

  g.summary();
}
}
