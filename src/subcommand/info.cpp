#include "./info.hpp"
#include "common.hpp"

namespace povu::subcommands::info {

void do_info(const core::config &app_config) {
  std::string fn_name = std::format("[povu::main::{}]", __func__);

  // -----
  // read the input gfa into a bidirected variation graph
  // -----

  bd::VG *g = povu::subcommands::common::get_vg(app_config);
  //bd::VG *g = io::from_gfa::to_bd(app_config.get_input_gfa().c_str(), app_config);

  g->summary();
  delete g;
}
}
