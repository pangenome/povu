#include "povu/prune.hpp"

#include <string>
#include <vector> // for vector

#include <log.h>		    // for log_info
#include <mto/from_gfa.hpp>	    // for to_bd
#include <mto/to_gfa.hpp>	    // for write_gfa
#include <quilt/app.hpp>	    // for config
#include <oza/graph/bidirected.hpp> // for bidirected
#include <oza/graph/pvst.hpp>	    // for pvst
#include <quilt/shim.hpp>	    // for format
#include <quilt/types.hpp>	    // for qt

namespace povu::subcommands::prune
{

void do_prune(const core::config &app_config)
{
	qt::u32 ll = app_config.verbosity();	      // ll for log level
	bd::VG *g = mto::from_gfa::to_bd(app_config); // read graph

	if (ll > 1)
		log_info("Finding components");

	std::vector<bd::VG *> components = bd::VG::componetize(*g);

	delete g;

	if (ll > 1)
		log_info("Found %ul components", components.size());

	std::string out_dir = app_config.get_output_dir();

	for (qt::u32 i{}; i < components.size(); ++i) {
		std::string fp =
			qs::format("{}/component_{}.gfa", out_dir, i + 1);
		mto::to_gfa::write_gfa(*components[i], fp);
	}

	return;
}

} // namespace povu::subcommands::prune
