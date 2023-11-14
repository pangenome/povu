#include <iostream>
#include <fstream>


#include "./core/core.hpp"
#include "./cli/cli.hpp"
#include "./io/io.hpp"
#include "./graph/bidirected.hpp"
#include "./graph/biedged.hpp"


int main(int argc, char *argv[]) {
	// set a higher value for tcmalloc warnings
	//setenv("TCMALLOC_LARGE_ALLOC_REPORT_THRESHOLD", "1000000000000000", 1);
  core::config app_config;
  cli::cli(argc, argv, app_config);

  app_config.dbg_print();

  bidirected::VariationGraph vg =
	io::from_gfa::to_vg(app_config.get_input_gfa().c_str());

  vg.dbg_print();

  biedged::BVariationGraph bg(vg);
  if (false) { bg.print_dot(); }
}
