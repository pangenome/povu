#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "./io/io.hpp"
#include "./pst/pst.hpp"
#include "./pvst/pvst.hpp"
#include "./graph/tree.hpp"
#include "./graph/digraph.hpp"
#include "./graph/u_graph.hpp"
#include "./graph/spanning_tree.hpp"
#include "./algorithms/cycle_equiv.hpp"
#include "./genomics/genomics.hpp"

// TODO: remove
#include "./example_graphs.hpp"
#include "core/core.hpp"


bool DEBUG = true;

void do_flubble(core::config& app_config) {
  
  // from a gfa file
  {
	
	digraph::DiGraph dg =
	  io::gfa_to_digraph( app_config.get_input_gfa().c_str());

	if (true) {
	  std::cout << "[domibubble::main] Successfully read GFA\n";
	  std::cout << "Graph properties:\n";
	}

	if (DEBUG) {
	  std::cout << "\n\n" << "Di-Graph" << "\n\n";
	  dg.print_dot();
	}

	u_graph::FlowGraph fg = u_graph::FlowGraph(dg);
	if (DEBUG) {
	  std::cout << "\n\n" << "Flow Graph" << "\n\n";
	  fg.print_dot();
	}

	spanning_tree::Tree t = fg.compute_spanning_tree();
	if (DEBUG) {
	  std::cout << "\n\n" << "Spanning tree" << "\n\n";
	  t.print_dot();
	}
	
	algorithms::cycle_equiv(t);
	if (DEBUG) {
	  std::cout << "\n\n" << "Updated Spanning tree" << "\n\n";
	  t.print_dot();
	}

	u_graph::FlowGraph afg = u_graph::FlowGraph(t);
	if (DEBUG) {
	  std::cout << "\n\n" << "Annotated Flow Graph" << "\n\n";
	  afg.print_dot();
	}

	std::vector<std::tuple< size_t , size_t, size_t>> v;
	std::vector<size_t> s;
	t.cycles_vector(v, s);

	tree::Tree pst_ =  pst::compute_pst(s);
	if (DEBUG) {
	  std::cout << "\n\n" << "PST" << "\n\n";
	  pst_.print_dot(true);
	}
	
	tree::Tree pvst_ =  pvst::compute_pvst(v);
	if (DEBUG) {
	  std::cout << "\n\n" << "PVST" << "\n\n";
	  pvst_.print_dot(true);
	}

	// call variants
	if (app_config.call_variants()) {
	  genomics::call_variants(pvst_, dg, app_config);
	}
	
	return;
  }

  {
  digraph::DiGraph g = g10();
  g = g20();
  
  if (DEBUG) {
    std::cout << "\n\n" << "Di-Graph" << "\n\n";
    g.print_dot();
  }

  if (DEBUG) {
    std::cout << "\n\n" << "Bi-edged Di-Graph" << "\n\n";
    //g.print_dot();
  }

  // a to_cfg() method
  //u_graph::FlowGraph u = g3();
  u_graph::FlowGraph u = u_graph::FlowGraph(g);
  if (DEBUG) {
    std::cout << "\n\n" << "Flow Graph" << "\n\n";
    u.print_dot();
  }

  spanning_tree::Tree t = u.compute_spanning_tree();
  if (DEBUG) {
    std::cout << "\n\n" << "Spanning tree" << "\n\n";
    t.print_dot();
  }

  algorithms::cycle_equiv(t);
  if (DEBUG) {
    std::cout << "\n\n" << "Updated Spanning tree" << "\n\n";
    t.print_dot();
  }

  std::vector<std::tuple< size_t , size_t, size_t>> v;
  std::vector<size_t> s;
  t.cycles_vector(v, s);

  //tree::Tree pst_ =  pvst::compute_pst(s);
  //pst_.print_dot(true);

  tree::Tree pvst_ =  pvst::compute_pvst(v);
  pvst_.print_dot(true);
  
  return;
  
  u_graph::FlowGraph afg = u_graph::FlowGraph(t);
  if (DEBUG) {
    std::cout << "\n\n" << "Annotated Flow Graph" << "\n\n";
    afg.print_dot();
  }

  spanning_tree::Tree sp2 = afg.compute_spanning_tree_two(g);
  if (DEBUG) {
    std::cout << "\n\n" << "Directed spanning tree Stack" << "\n\n";
    sp2.print_dot();

  }

  sp2.compute_edge_stack();

  return;
  
  std::vector<u_graph::Edge> edge_stack = afg.compute_edge_stack();
  if (DEBUG) {
    std::cout << "\n\n" << "Edge Stack" << "\n\n";
    for (auto e : edge_stack) {
      std::cout << "" << e.left()  << "-- (" << e.get_weight()  << ") --" << e.right() << "\n";
    }
  }

  tree::Tree pst =  afg.construct_pst(edge_stack);
  if (DEBUG) {
    std::cout << "\n\n" << "PST" << "\n\n";
    pst.print_dot(true);
  }

  tree::Tree pvst =  afg.construct_pvst(edge_stack);
  if (DEBUG) {
    std::cout << "\n\n" << "PVST" << "\n\n";
    pvst.print_dot(true);
  }

  return;
  }
}

int main() {

  core::config app_config;
	
  std::string filename = "./test_data/t3.gfa";
  app_config.set_input_gfa(filename.c_str());
	
  app_config.add_reference_path("hap1");

  app_config.set_call_variants(true);
  
  do_flubble(app_config);

  return 0;
}
