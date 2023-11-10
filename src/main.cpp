#include <iostream>
#include <fstream>


#include "./core/core.hpp"
#include "./cli/cli.hpp"
#include "./io/io.hpp"
#include "./graph/bidirected.hpp"


int main(int argc, char *argv[]) {
    // set a higher value for tcmalloc warnings
    //setenv("TCMALLOC_LARGE_ALLOC_REPORT_THRESHOLD", "1000000000000000", 1);
  core::config app_config;
  cli::cli(argc, argv, app_config);

  std::cout << "input_gfa: " << app_config.get_input_gfa() << std::endl;
  std::cout << "l:" << app_config.get_reference_paths().size();
  for (auto &&path : app_config.get_reference_paths()) {
	std::cout << "path: " << path << std::endl;

	
  }


  bidirected::VariationGraph vg =
  io::gfa_to_vg(app_config.get_input_gfa().c_str());

  

}
