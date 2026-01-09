#include "./prune.hpp"

#include <string>
#include <vector> // for vector

#include "povu/common/app.hpp"	     // for config
#include "povu/common/compat.hpp"    // for pv_cmp, format
#include "povu/graph/bidirected.hpp" // for bidirected
#include "povu/graph/pvst.hpp"	     // for pvst
#include "povu/io/from_gfa.hpp"	     // for to_bd
#include "povu/io/to_gfa.hpp"	     // for write_gfa

namespace povu::subcommands::prune
{

void do_prune(const core::config &app_config)
{
	pt::u32 ll = app_config.verbosity();		   // ll for log level
	bd::VG *g = povu::io::from_gfa::to_bd(app_config); // read graph

	if (ll > 1)
		INFO("Finding components");

	std::vector<bd::VG *> components = bd::VG::componetize(*g);

	delete g;

	if (ll > 1)
		INFO("Found {} components", components.size());

	std::string out_dir = app_config.get_output_dir();

	for (pt::u32 i{}; i < components.size(); ++i) {
		std::string fp =
			pv_cmp::format("{}/component_{}.gfa", out_dir, i + 1);
		povu::io::to_gfa::write_gfa(*components[i], fp);
	}

	return;
}

} // namespace povu::subcommands::prune
