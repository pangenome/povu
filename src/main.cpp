#include <algorithm>
#include <cassert>
#include <cstddef>
#include <format>
#include <iostream>
#include <string>
#include <stack>
#include <tuple>
#include <utility>


#include "./algorithms/cycle_equiv.hpp"
#include "./cli/cli.hpp"
#include "./core/core.hpp"
#include "./genomics/genomics.hpp"
#include "./graph/bidirected.hpp"
#include "./graph/biedged.hpp"
#include "./graph/spanning_tree.hpp"



#include "./io/io.hpp"

#include "./common/common.hpp"
#include "core/constants.hpp"


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

  if (app_config.verbosity() > 1) {
    if (app_config.verbosity() > 2)  { std::cerr << fn_name << " Finished reading graph:\n"; }
    // no need to print this if we have multiple components
    // vg.dbg_print();
  }

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

spanning_tree::Tree compute_spanning_tree(biedged::BVariationGraph g) {
  spanning_tree::Tree t = spanning_tree::Tree();

  std::stack<std::tuple<std::size_t, std::size_t, bidirected::color>> s; // parent and v idx
  std::set<std::size_t> visited;
  std::size_t v_idx {}; // set start node to 0
  s.push({core::constants::INVALID_ID, v_idx, bidirected::color::gray});
  visited.insert(v_idx);
  std::size_t counter{}; // dfs pre visit counter

  std::map<std::size_t, std::size_t> vtx_to_dfs_num;
  std::set<std::size_t> in_tree; // graph vertices in the tree
  std::set<std::pair<std::size_t, std::size_t>> back_edges; // avoids duplicate back edges

  auto is_parent = [&](std::size_t p, std::size_t c) -> bool {
    return t.get_children(vtx_to_dfs_num[p]).count(vtx_to_dfs_num[c]);
  };

  auto is_parent_child = [&](std::size_t n, std::size_t v_idx) -> bool {
    return is_parent(n, v_idx) || is_parent(v_idx, n);
  };

  while (!s.empty()) {
    auto [p_idx, v_idx, c] = s.top();

    biedged::Vertex const& v = g.get_vertex(v_idx);

    if (in_tree.find(v_idx) == in_tree.end()) {
      t.add_vertex(spanning_tree::Vertex{v_idx, counter, v.get_handle(), v.get_type()});
      in_tree.insert(v_idx);
      vtx_to_dfs_num[v_idx] = counter;

      if (p_idx != core::constants::INVALID_ID) {
        t.add_tree_edge(vtx_to_dfs_num[p_idx], counter, core::constants::UNDEFINED_SIZE_T, c);
      }

      counter++;
    }

    bool explored {true};

    //std::set<std::pair<graph_types::color, std::size_t>> neighbors =
    for (auto [c, n] : g.get_neighbours(v_idx)) {
      if (visited.find(n) == visited.end()) {
        s.push({v_idx, n, c});
        visited.insert(n);
        explored = false;
        break;
      }
      else if (!is_parent_child(n, v_idx) && back_edges.find({n, v_idx}) == back_edges.end())  {
        t.add_be(vtx_to_dfs_num[v_idx], vtx_to_dfs_num[n],core::constants::UNDEFINED_SIZE_T, false, c);
        back_edges.insert({v_idx, n});
      }
    }

    if (explored) { s.pop(); }
  }

  // TODO: wrap in DEBUG pragma
  // check the correctness of the tree
  {

    if (t.size() != g.size()) {
      std::cerr << "size mismatch " << t.size() << " " << g.size() << "\n";
    }

    // TODO: compare edge count

  for (std::size_t v{}; v < t.size() ; v++) {
    for (auto c : t.get_children(v)) {
      if (t.get_vertex(v).dfs_num() >= t.get_vertex(c).dfs_num()) {
        std::cerr << "te weird " << v << " " << t.get_vertex(v).dfs_num() << " "
                  << c << " " << t.get_vertex(c).dfs_num() << "\n";
      }
    }


    for (std::size_t tgt : t.get_obe(v)) {
      if (t.get_vertex(v).dfs_num() <= t.get_vertex(tgt).dfs_num()) {
        std::cerr << "obe weird " << v << " " << t.get_vertex(v).dfs_num() << " "
                  << tgt << " " << t.get_vertex(tgt).dfs_num() << "\n";
      }
    }

    for (std::size_t src : t.get_ibe(v)) {
      if (t.get_vertex(v).dfs_num() >= t.get_vertex(src).dfs_num()) {
         std::cerr << "ibe weird" << v << " " << t.get_vertex(v).dfs_num() << " "
                   << src << " " << t.get_vertex(src).dfs_num() << "\n";
      }
    }
  }
}
  return t;
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
  //spanning_tree::Tree st =    bg.compute_spanning_tree();

  spanning_tree::Tree st = compute_spanning_tree(bg);


  if (app_config.print_dot() && app_config.verbosity() > 4) { std::cout << "\n\n" << "Spanning Tree" << "\n\n";
    st.print_dot();
  }





  if (app_config.verbosity() > 2) { std::cerr << fn_name << " Finding Cycle Equivalent Classes\n"; }
  algorithms::cycle_equiv(st);

  return;


  if (app_config.print_dot() && app_config.verbosity() > 4) { std::cout << "\n\n" << "Updated Spanning Tree" << "\n\n";
    st.print_dot();
  }

  return;

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

  for (std::size_t i{8}; i < components.size(); i++) {
    if (app_config.verbosity() > 2) {
      std::cerr << std::format("{} Handling component: {}\n", fn_name, i+1);
    }

    if (components[i].size() < 2) {
      std::cerr << std::format("{} Skipping component because size {} is too small\n", fn_name, components[i].size());
      continue;
    }

    if (app_config.verbosity() > 3) { components[i].dbg_print(); }


    // if ( i==8) { continue; }

    compute_sese_regions(components[i], app_config);

    return 0;
  }

  return 0;
}
