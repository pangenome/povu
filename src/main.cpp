#include <string>
#include <format>

#include "./algorithms/cycle_equiv.hpp"
#include "./cli/cli.hpp"
#include "./core/core.hpp"
#include "./genomics/genomics.hpp"
#include "./graph/bidirected.hpp"
#include "./graph/biedged.hpp"
#include "./graph/spanning_tree.hpp"
#include "./graph/tree.hpp"
#include "./graph/u_graph.hpp"
#include "./io/io.hpp"
#include "./pvst/pvst.hpp"


int main(int argc, char *argv[]) {
  std::string fn_name = std::format("[povu::main::{}]", __func__);

  core::config app_config;
  cli::cli(argc, argv, app_config);

  if (app_config.verbosity()) { app_config.dbg_print(); }

  if (app_config.verbosity() > 2)  { std::cerr << fn_name << " Reading graph\n"; }

  // read the input gfa into a bidirected variation graph
  bidirected::VariationGraph vg = io::from_gfa::to_vg(app_config.get_input_gfa().c_str(), app_config);
  if (app_config.verbosity() > 1) { vg.dbg_print(); }

  if (app_config.verbosity() > 2) { std::cerr << fn_name << " Bi-edging" << "\n"; }
  // convert the bidirected variation graph into a biedged variation graph
  biedged::BVariationGraph bg(vg);
  if (app_config.verbosity() > 4) { bg.print_dot(); }
  bg.componetize();
  if (app_config.verbosity() > 4) { bg.print_dot(); }

  // compute the spanning tree of the biedged variation graph
  if (app_config.verbosity() > 2) { std::cerr << fn_name << " Computing spanning tree\n"; }
  spanning_tree::Tree st = bg.compute_spanning_tree();
  if (app_config.verbosity() > 4) { std::cout << "\n\n" << "Spanning Tree" << "\n\n"; st.print_dot(); }

  // compute the cycle equivalence classes of the spanning tree
  if (app_config.verbosity() > 2) { std::cerr << fn_name << " Finding cycle equivalent classes\n"; }
  algorithms::cycle_equiv(st);
  if (app_config.verbosity() > 4) {
    std::cout << "\n\n" << "Annotated Spanning Tree" << "\n\n";
    st.print_dot();
  }

  std::vector<core::eq_n_id_t> v = st.compute_edge_stack2();

  u_graph::FlowGraph afg = u_graph::FlowGraph(st);
  if (app_config.verbosity() > 2) {
    std::cout << "\n\n" << "Annotated Flow Graph" << "\n\n";
    afg.print_dot();
  }

  if (app_config.verbosity() > 2) { std::cerr << fn_name << " Computing PVST\n"; }
  tree::Tree t = pvst::compute_pvst(v, app_config);
  if (app_config.verbosity() > 4)  { std::cout << "\n\n" << "PVST" << "\n\n"; t.print_dot(true); }

  //return 0;

  if (app_config.verbosity() > 2)  { std::cerr << fn_name << " Calling variants\n"; }
  genomics::call_variants(t, vg, app_config);


  //return 0;


  return 0;
}
