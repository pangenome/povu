#include <iostream>

#include "./io/io.hpp"
#include "./algorithms/cycle_equiv.hpp"
//#include "./vst/vst.hpp"
//#include "./vst/pst.hpp"
#include "./pvst/pvst.hpp"
#include "./pst/pst.hpp"
#include "./graph/tree.hpp"
#include "./graph/digraph.hpp"
#include "./graph/u_graph.hpp"
#include "./graph/spanning_tree.hpp"

#include "./example_graphs.hpp"


bool DEBUG = true;

int main() {

  // from a gfa file
  {
	digraph::DiGraph dg;
	io::gfa_to_digraph("/home/sluggie/src/phd/domibubble-cpp/deps/gfakluge/data/test.gfa", &dg);
	dg.print_dot();

	u_graph::FlowGraph fg = u_graph::FlowGraph(dg);
	if (DEBUG) {
	  std::cout << "\n\n" << "Flow Graph" << "\n\n";
	  fg.print_dot();
	}
  
	return 0;
  }
  
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
  
  return 0;
  
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

  return 0;
  
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

  return  0;
}
