#include <algorithm>
#include <cassert>
#include <cstddef>
#include <format>
#include <iostream>
#include <string>
#include <utility>

#include "./algorithms/cycle_equiv.hpp"
#include "./cli/cli.hpp"
#include "./common/common.hpp"
#include "./cli/app.hpp"
#include "./graph/bidirected.hpp"
#include "./graph/biedged.hpp"
#include "./io/io.hpp"



/**
 * @brief reads the input gfa into a bidirected variation graph, and returns the components
 *
 * Given a bidirected::VariationGraph return each component as a map of id_t to component
 * each component is a bidirected::component which contains a bidirected::VariationGraph
 *
 * @param app_config
 * @return std::map<id_t, bidirected::component>
 */
std::vector<bidirected::VariationGraph> read_and_componetize(const core::config& app_config) {
  std::string fn_name = std::format("[povu::main::{}]", __func__);

  if (app_config.verbosity() > 2)  { std::cerr << fn_name << " Reading graph\n"; }

  // read the input gfa into a bidirected variation graph
  bidirected::VariationGraph vg =
    io::from_gfa::to_vg(app_config.get_input_gfa().c_str(), app_config);

  /*
    no need to print this if we have multiple components
    if (app_config.verbosity() > 1) {
    if (app_config.verbosity() > 2)  { std::cerr << fn_name << " Finished reading graph:\n"; }
    // vg.dbg_print();
    }
   */

  if (false) { // validate the haplotype paths
    if (app_config.verbosity() > 2)  { std::cerr << fn_name << " Validating paths\n"; }
    if (vg.validate_haplotype_paths()) {
      std::cerr << fn_name << " Haplotype paths are valid" << "\n";
    }
  }

  return bidirected::componetize(vg, app_config);
}

void update_bd_eq_classes(const biedged::BVariationGraph &be_g, bidirected::VariationGraph &bd_g) {
  std::string fn_name = std::format("[povu::main::{}]", __func__);

  std::vector<std::size_t> d_vs = be_g.get_dummy_vertices();
  auto bidirected_idx = [&](std::size_t x) -> std::size_t {
    if (!d_vs.empty()) { --x; } // we added 1 because we added a dummy start node before bi-edging
    return x % 2 == 0 ? ((x + 2) / 2) - 1 : ((x + 1) / 2) - 1;
  };

  const std::vector<bidirected::Edge>& bidirected_edges = bd_g.get_all_edges();
  // a map of vertex pair to edge index in the bidirected graph
  std::map<std::pair<std::size_t, std::size_t>, std::size_t> vertex_pair_to_edge_idx;
  for (std::size_t i = 0; i < bidirected_edges.size(); ++i) {
    const bidirected::Edge& e = bidirected_edges[i];
    std::size_t v1 = e.get_v1_idx();
    std::size_t v2 = e.get_v2_idx();
    v1 = std::stoull(bd_g.get_vertex(v1).get_name());
    v2 = std::stoull(bd_g.get_vertex(v2).get_name());
    vertex_pair_to_edge_idx[std::make_pair(v1, v2)] = i;
    vertex_pair_to_edge_idx[std::make_pair(v2, v1)] = i;
  }

  // inline a function to take two vertices and return whether they are in d_vs or not
  auto is_dummy = [&](std::size_t v1, std::size_t v2) -> bool {
    return std::find(d_vs.begin(), d_vs.end(), v1) != d_vs.end() ||
           std::find(d_vs.begin(), d_vs.end(), v2) != d_vs.end();
  };

  const std::vector<biedged::Edge>& biedged_edges = be_g.get_all_edges();

  for (const biedged::Edge& e : biedged_edges) {
    std::size_t v1 = e.get_v1_idx();
    std::size_t v2 = e.get_v2_idx();

    // skip the dummy vertices
    if (is_dummy(v1, v2)) { continue; }

    graph_types::color c = e.get_color();

    if (c == graph_types::color::black) {

      std::size_t bd_v1_idx = bidirected_idx(v1);
      std::size_t bd_v2_idx = bidirected_idx(v2);

      assert(bd_v1_idx == bd_v2_idx);

      bd_g.get_vertex_mut(bd_v1_idx).set_eq_class(e.get_eq_class());
      bd_g.get_vertex_mut(bd_v2_idx).set_eq_class(e.get_eq_class());
    }
    else if (c == graph_types::color::gray) {
      v1 = std::stoull(be_g.get_vertex(v1).get_handle());
      v2 = std::stoull(be_g.get_vertex(v2).get_handle());

      std::size_t bd_e_idx = vertex_pair_to_edge_idx[std::make_pair(v1, v2)];
      bd_g.get_edge_mut(bd_e_idx).set_eq_class(e.get_eq_class());
    }
  }
}


/**
 * @brief takes a variation graph which should be in a single component and computes its SESE regions
 *
 * @param app_config the configuration
 * @param vg the variation graph
 */
void compute_sese_regions(bidirected::VariationGraph &vg, core::config& app_config) {
  std::string fn_name = std::format("[povu::main::{}]", __func__);

  if (app_config.print_dot() && app_config.verbosity() > 4) { std::cout << "\n\n" << "Variation Graph (unsorted)" << "\n\n";
    vg.print_dot();
  }

  // convert the bidirected variation graph into a biedged variation graph
  if (app_config.verbosity() > 2) { std::cerr << fn_name << " Bi-edging" << "\n"; }
  biedged::BVariationGraph bg(vg); // will add dummy vertices
  if (app_config.print_dot() && app_config.verbosity() > 4 ) { std::cout << "\n\n" << "Biedged" << "\n\n";
    bg.print_dot();
  }

  // compute the spanning tree of the biedged variation graph
  if (app_config.verbosity() > 2) { std::cerr << fn_name << " Generating spanning tree\n"; }
  spanning_tree::Tree st = bg.compute_spanning_tree();

  if (app_config.print_dot() && app_config.verbosity() > 4) { std::cout << "\n\n" << "Spanning Tree" << "\n\n";
    st.print_dot();
  }

  if (app_config.verbosity() > 2) { std::cerr << fn_name << " Computing Cycle Equivalence\n"; }
  algorithms::cycle_equiv(st);

  if (app_config.print_dot() && app_config.verbosity() > 4) { std::cout << "\n\n" << "Updated Spanning Tree" << "\n\n";
    st.print_dot();
  }

  if (app_config.verbosity() > 2) { std::cerr << fn_name << " Finding SESE regions\n"; }
  algorithms::find_seses(st);

  return;
  /*
  bg.update_eq_classes(st);
  if (app_config.print_dot() && app_config.verbosity() > 4) { std::cout << "\n\n" << "Updated Biedged" << "\n\n";
    bg.print_dot();
  }

  update_bd_eq_classes(bg, vg);

  if (app_config.print_dot() && app_config.verbosity() > 4) { std::cout << "\n\n" << "Variation Graph (with classes, unsorted)" << "\n\n";
    vg.print_dot();
  }

  if (app_config.verbosity() > 2) { std::cerr << fn_name << " Computing eq class stack\n"; }
  std::vector<std::size_t> v = algorithms::compute_eq_class_stack(st);

  if (app_config.verbosity() > 2) { std::cerr << fn_name << " Finding SNPs\n"; }
  std::vector<std::tuple<std::size_t, graph_types::VertexType, std::size_t, graph_types::VertexType>> seses =
    graph_operations::foo(st, v, app_config);

  if (app_config.verbosity() > 2) { std::cerr << fn_name << " SESE regions\n"; }
  for (const auto& [v1, t1, v2, t2] : seses) {

    std::string n1 {vg.get_vertex(v1).get_name()};
    std::string n2 {vg.get_vertex(v2).get_name()};
    //std::cout << "(" << n1 << ", " << t1 << ") ~> (" << n2 << ", " << t2 << ")\n";
  }
  */
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

  std::vector<bidirected::VariationGraph> components = read_and_componetize(app_config);

  if (app_config.verbosity() > 2)  {
    std::cerr << std::format("{} Number of components: {}\n", fn_name, components.size());
  }

  for (std::size_t i{}; i < components.size(); i++) {
    if (app_config.verbosity() > 2) {
      std::cerr << std::format("{} Handling component: {}\n", fn_name, i+1);
    }

    if (components[i].size() < 2) {
      std::cerr << std::format("{} Skipping component because size {} is too small\n", fn_name, components[i].size());
      continue;
    }

    if (app_config.verbosity() > 3) { components[i].dbg_print(); }

    if (i==6) {
      std::cerr << "Skipping component 6\n";
      continue;
    }

    compute_sese_regions(components[i], app_config);
  }

  return 0;
}
