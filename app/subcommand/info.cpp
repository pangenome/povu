#include "./info.hpp"

#include <iostream>                   // for basic_ostream, basic_ios, cerr
#include <string>                     // for basic_string, operator<<, string
#include <vector>                     // for vector

#include "povu/graph/bidirected.hpp"  // for VG
#include "povu/io/from_gfa.hpp"       // for to_bd
#include "fmt/core.h"                 // for format
#include "povu/common/compat.hpp"     // for format, pv_cmp
#include "povu/common/core.hpp"       // for idx_t, pt

namespace povu::subcommands::info
{
namespace bd = povu::bidirected;

void do_info(const core::config &app_config)
{
	std::string fn_name = pv_cmp::format("[povu::main::{}]", __func__);

	// -----
	// read the input gfa into a bidirected variation graph
	// -----

	bd::VG *g = povu::io::from_gfa::to_bd(app_config);

	std::vector<bd::VG *> components = bd::VG::componetize(*g);

	delete g;

	std::cerr << pv_cmp::format("{} Component count {}\n", fn_name,
				    components.size());

	for (pt::idx_t i{}; i < components.size(); ++i) {
		bd::VG *c = components[i];
		c->summary(app_config.print_tips());
		delete c;
		components[i] = nullptr;
	}
	components.clear();

	return;
}
} // namespace povu::subcommands::info
