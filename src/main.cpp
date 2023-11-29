#include "./core/core.hpp"
#include "./cli/cli.hpp"
#include "./io/io.hpp"
#include "./graph/bidirected.hpp"
#include "./graph/biedged.hpp"
#include "./graph/spanning_tree.hpp"
#include "./algorithms/cycle_equiv.hpp"
#include "./pvst/pvst.hpp"

#include "./graph/u_graph.hpp"
#include "graph/tree.hpp"

int main(int argc, char *argv[]) {
  core::config app_config;
  cli::cli(argc, argv, app_config);

  if (app_config.verbosity()) { app_config.dbg_print(); }

  // read the input gfa into a bidirected variation graph
  bidirected::VariationGraph vg = io::from_gfa::to_vg(app_config.get_input_gfa().c_str());
  if (app_config.verbosity() > 1) { vg.dbg_print(); }

  // convert the bidirected variation graph into a biedged variation graph
  biedged::BVariationGraph bg(vg);
  if (app_config.verbosity() > 2) { bg.print_dot(); }
  bg.componetize();
  if (app_config.verbosity() > 2) { bg.print_dot(); }

  // compute the spanning tree of the biedged variation graph
  spanning_tree::Tree st = bg.compute_spanning_tree();
  if (app_config.verbosity() > 2) {
	std::cout << "\n\n" << "Spanning Tree" << "\n\n";
	st.print_dot();
  }

  // compute the cycle equivalence classes of the spanning tree
  algorithms::cycle_equiv(st);
  if (app_config.verbosity() > 2) {

	st.print_dot();
  }

  std::vector<std::pair<std::size_t, std::size_t>> v = st.compute_edge_stack2();

  tree::Tree t = pvst::compute_pvst(v);
  if (app_config.verbosity() > 2)  {
	std::cout << "\n\n" << "PVST" << "\n\n";
	t.print_dot(true);	
  }


  return 0;
  
  u_graph::FlowGraph afg = u_graph::FlowGraph(st);
  if (app_config.verbosity() > 2) { 
    std::cout << "\n\n" << "Annotated Flow Graph" << "\n\n";
    afg.print_dot();
  }
  

  return 0;
}
