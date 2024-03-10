#include <algorithm>
#include <cstddef>
#include <format>
#include <string>

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


/**
  * @brief reads the input gfa into a bidirected variation graph, and returns the components
  *
  * Given a bidirected::VariationGraph return each component as a map of id_t to component
  * each component is a bidirected::component which contains a bidirected::VariationGraph
  *
  * @param app_config
  * @return std::map<id_t, bidirected::component>
 *
 */
std::map<id_t, bidirected::component> read_and_componetize(const core::config& app_config) {
  std::string fn_name = std::format("[povu::main::{}]", __func__);

  if (app_config.verbosity() > 2)  { std::cerr << fn_name << " Reading graph\n"; }

  // read the input gfa into a bidirected variation graph
  bidirected::VariationGraph vg =
    io::from_gfa::to_vg(app_config.get_input_gfa().c_str(), app_config);

  if (app_config.verbosity() > 1) {
    vg.dbg_print();
  }

  if (app_config.verbosity() > 4) { std::cout << "\n\n" << "Variation Graph (unsorted)" << "\n\n";
    vg.print_dot();
  }

  std::map<id_t, bidirected::component> components = vg.count_components(app_config);
  return components;
}

/**
  * @brief takes a variation graph which should be in a single component and returns the
  *        calls variants on it
  *
  * @param app_config the configuration
  * @param vg the variation graph
  */
void bar(core::config& app_config, bidirected::VariationGraph vg) {
  std::string fn_name = std::format("[povu::main::{}]", __func__);

  if (app_config.verbosity() > 1) {
    vg.dbg_print();
  }

  if (app_config.verbosity() > 2)  { std::cerr << fn_name << " Sorting graph\n"; }


  if (app_config.sort()) {
    vg.sort();
  }



  if (app_config.verbosity() > 1) {
    vg.dbg_print();
  }

  if (app_config.verbosity() > 4) { std::cout << "\n\n" << "Variation Graph (sorted)" << "\n\n";
    vg.print_dot();
  }

  if (app_config.verbosity() > 2) { std::cerr << fn_name << " Bi-edging" << "\n"; }
  // convert the bidirected variation graph into a biedged variation graph
  biedged::BVariationGraph bg(vg);
  if (app_config.verbosity() > 4) {
    std::cout << "\n\n" << "Biedged" << "\n\n";
    bg.print_dot();
  }
  bg.componetize();
  if (app_config.verbosity() > 4) { std::cout << "\n\n" << "Componetized biedged" << "\n\n";
    bg.print_dot();
  }

  // compute the spanning tree of the biedged variation graph
  if (app_config.verbosity() > 2) { std::cerr << fn_name << " Computing spanning tree\n"; }
  spanning_tree::Tree st = bg.compute_spanning_tree();
  if (app_config.verbosity() > 4) { std::cout << "\n\n" << "Spanning Tree" << "\n\n";
    st.print_dot();
  }

  // compute the cycle equivalence classes of the spanning tree
  if (app_config.verbosity() > 2) { std::cerr << fn_name << " Finding cycle equivalent classes\n"; }
  algorithms::cycle_equiv(st);
  if (app_config.verbosity() > 4) { std::cout << "\n\n" << "Annotated Spanning Tree" << "\n\n";
    st.print_dot();
  }

  std::vector<core::eq_n_id_t> v = st.compute_edge_stack2();

  u_graph::FlowGraph afg = u_graph::FlowGraph(st);
  if (app_config.verbosity() > 2) { std::cout << "\n\n" << "Annotated Flow Graph" << "\n\n";
    afg.print_dot();
  }

  if (app_config.verbosity() > 2) { std::cerr << fn_name << " Computing PVST\n"; }
  tree::Tree t = pvst::compute_pvst(v, app_config);
  if (app_config.verbosity() > 2) { std::cout << "\n\n" << "PVST" << "\n\n";
    t.print_dot(true);
  }

  if (app_config.get_pvst_path().has_value()) {
    pvst::to_text(t, app_config.get_pvst_path().value());
  }

  if (app_config.verbosity() > 2)  { std::cerr << fn_name << " Calling variants\n"; }
  genomics::call_variants(t, vg, app_config);

  return;
}

/**
 * @brief main function
 *
 * @param argc
 * @param argv
 * @return int
 */
int main(int argc, char *argv[]) {
  std::string fn_name = std::format("[povu::main::{}]", __func__);

  core::config app_config;
  cli::cli(argc, argv, app_config);

  if (app_config.verbosity()) { app_config.dbg_print(); }

  std::map<id_t, bidirected::component> components = read_and_componetize(app_config);

  if (app_config.verbosity() > 2)  {
    std::cerr << std::format("{} Number of components: {}\n", fn_name, components.size());
  }

  for (auto const &[_, v] : components) {
    bar(app_config, v.vg);
    break;
  }

  return 0;
}
