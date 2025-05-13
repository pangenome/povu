#include "./info.hpp"
#include "common.hpp"

namespace povu::subcommands::info {

void do_info(const core::config &app_config) {
  std::string fn_name = std::format("[povu::main::{}]", __func__);

  // -----
  // read the input gfa into a bidirected variation graph
  // -----

  bd::VG *g = povu::subcommands::common::get_vg(app_config);

  std::vector<bd::VG *> components = bd::componetize(*g);

  delete g;

  std::cerr << std::format("{} Component count {}\n", fn_name, components.size());

  for (pt::idx_t i{}; i < components.size(); ++i) {
    bd::VG *c = components[i];
    c->summary(app_config.print_tips());
    delete c;
    components[i] = nullptr;
  }
  components.clear();

  return;
}
}
