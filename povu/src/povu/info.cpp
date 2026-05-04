#include "povu/info.hpp"

#include <iostream> // for basic_ostream, basic_ios, cerr
#include <string>   // for basic_string, operator<<, string
#include <vector>   // for vector

#include <mto/from_gfa.hpp>	    // for to_bd
#include <oza/graph/bidirected.hpp> // for VG
#include <quilt/shim.hpp>	    // for format
#include <quilt/types.hpp>	    // for qt

namespace povu::subcommands::info
{

void do_info(const core::config &app_config)
{
	std::string fn_name = qs::format("[povu::main::{}]", __func__);

	// -----
	// read the input gfa into a bidirected variation graph
	// -----

	bd::VG *g = mto::from_gfa::to_bd(app_config);

	std::vector<bd::VG *> components = bd::VG::componetize(*g);

	delete g;

	std::cerr << qs::format("{} Component count {}\n", fn_name,
				components.size());

	for (qt::idx_t i{}; i < components.size(); ++i) {
		bd::VG *c = components[i];
		c->summary(app_config.print_tips());
		delete c;
		components[i] = nullptr;
	}
	components.clear();

	return;
}
} // namespace povu::subcommands::info
