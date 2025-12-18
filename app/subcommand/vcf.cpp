#include "./vcf.hpp"

#include <string>
#include <vector> // for vector

#include "povu/common/app.hpp"	     // for config
#include "povu/common/compat.hpp"    // for pv_cmp, format
#include "povu/graph/bidirected.hpp" // for bidirected
#include "povu/graph/pvst.hpp"	     // for pvst
#include "povu/io/from_gfa.hpp"	     // for to_bd
#include "povu/io/to_gfa.hpp"	     // for write_gfa

namespace povu::subcommands::vcf
{

void do_vcf(const core::config &app_config)
{
	pt::u32 ll = app_config.verbosity();		   // ll for log level
	bd::VG *g = povu::io::from_gfa::to_bd(app_config); // read graph

	delete g;

	return;
}

} // namespace povu::subcommands::vcf
