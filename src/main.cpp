#include <iostream>

#include "./graph/digraph.hpp"
#include "./graph/graph.hpp"
#include "./graph/spanning_tree.hpp"
#include "./vst/vst.hpp"
#include "./vst/pst.hpp"

u_graph::CFG g1() {
  u_graph::CFG g;

  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(2, 3);
  g.add_edge(3, 1);

  g.set_start_node(0);
  g.set_stop_node(3);

  return g;
};

// big
u_graph::CFG g2() {
  u_graph::CFG g;
  g = u_graph::CFG();

  g.add_edge(0, 3);
  g.add_edge(0, 1);
  g.add_edge(3, 2);
  g.add_edge(1, 2);

  g.add_edge(2, 9);
  g.add_edge(2, 8);

  g.add_edge(1, 4);
  g.add_edge(1, 6);
  g.add_edge(1, 7);

  g.add_edge(5, 4);
  g.add_edge(7, 4);

  g.add_edge(4, 6);
  g.add_edge(7, 6);

  g.set_start_node(0);
  g.set_stop_node(5);
  g.set_stop_node(8);
  g.set_stop_node(9);

  return g;
};

// diamond
u_graph::CFG g3() {
  u_graph::CFG g;
  g = u_graph::CFG();

  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(1, 3);
  g.add_edge(2, 4);
  g.add_edge(3, 4);
  g.add_edge(4, 5);

  g.set_start_node(0);
  g.set_stop_node(5);

  return g;
}

// overleaf fig 1
// missed by bubble finder
u_graph::CFG g4() {
  u_graph::CFG g;
  g = u_graph::CFG();

  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(1, 3);
  g.add_edge(2, 4);
  g.add_edge(3, 4);
  g.add_edge(1, 5);
  g.add_edge(4, 5);

  g.set_start_node(0);
  g.set_stop_node(5);

  return g;
}

// fig 1 a paper
// diamond
u_graph::CFG g5() {
  u_graph::CFG g;
  g = u_graph::CFG();

  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(1, 3);
  g.add_edge(2, 4);
  g.add_edge(3, 4);
  g.add_edge(4, 5);
  g.add_edge(3, 5);
  g.add_edge(6, 5);

  g.set_start_node(0);
  g.set_stop_node(6);

  return g;
}

// overlapping
u_graph::CFG g6() {
  u_graph::CFG g;

  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(2, 3);
  g.add_edge(3, 4);
  g.add_edge(4, 2);
  g.add_edge(4, 5);
  g.add_edge(1, 5);
  g.add_edge(5, 6);

  g.set_start_node(0);
  g.set_stop_node(6);

  return g;
};

digraph::DiGraph g7() {
  digraph::DiGraph g(7);

  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(2, 3);
  g.add_edge(3, 4);
  g.add_edge(4, 2);
  g.add_edge(4, 5);
  g.add_edge(1, 5);
  g.add_edge(5, 6);

  g.add_edge(5, 7);
  g.add_edge(6, 7);

  g.add_edge(7, 8);

  g.add_start_node(0);
  g.add_stop_node(6);

  return g;

}

digraph::DiGraph g8() {
  digraph::DiGraph g(5);

  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(0, 2);
  g.add_edge(2, 3);
  g.add_edge(2, 4);

  g.add_start_node(0);
  g.add_stop_node(3);
  g.add_stop_node(4);

  return g;
}

// from ekg
digraph::DiGraph g9() {
  digraph::DiGraph g(16);

  g.add_edge(0, 1);

  g.add_edge(1, 3);
  g.add_edge(3, 5);

  g.add_edge(5, 7);
  g.add_edge(5, 9);
  g.add_edge(7, 11);
  g.add_edge(9, 11);
  g.add_edge(11, 13);
  g.add_edge(13, 1);

  g.add_edge(0, 2);
  g.add_edge(2, 4);
  g.add_edge(4, 6);
  g.add_edge(4, 8);
  g.add_edge(6, 10);
  g.add_edge(10, 8);
  g.add_edge(10, 12);
  g.add_edge(8, 12);
  g.add_edge(12, 14);
  g.add_edge(14, 15);
  g.add_edge(15, 0);
  g.add_edge(13, 15);

  g.add_start_node(0);
  g.add_stop_node(15);

  return g;
}

// non planar 1
digraph::DiGraph g10() {
  digraph::DiGraph g(5);

  g.add_edge(0, 1);

  g.add_edge(1, 3);
  g.add_edge(1, 2);

  g.add_edge(3, 4);
  g.add_edge(2, 4);
  g.add_edge(0, 4);
  g.add_edge(3, 5);
    g.add_edge(2, 5);

  g.add_start_node(0);
  g.add_stop_node(5);

  return g;
}

// non planar 2
digraph::DiGraph g11() {
  digraph::DiGraph g(5);

  g.add_edge(0, 1);
  g.add_edge(0, 3);
  g.add_edge(1, 3);
  g.add_edge(1, 2);
  g.add_edge(2, 3);
  g.add_edge(2, 4);

  g.add_start_node(0);
  g.add_stop_node(4);

  return g;
}

digraph::DiGraph g12() {
  digraph::DiGraph g;

  g.add_edge(0, 1);
    g.add_edge(0, 4);
  g.add_edge(1, 2);
  g.add_edge(1, 3);

  g.add_edge(2, 3);
  g.add_edge(2, 4);

  g.add_start_node(0);
  g.add_stop_node(4);

  return g;
};


const bool DEBUG = true;

int main() {

  digraph::DiGraph g = g10();

  g.print_dot();
  g.biedge();
  g.print_dot();

  // a to_cfg() method
  u_graph::CFG u = u_graph::CFG(g);
  u.print_dot();


  spanning_tree::Tree t = u.compute_spanning_tree();
  if (DEBUG) {
    std::cout << "\n\n" << "Spanning tree" << "\n\n";
    t.print_dot();
  }

  /*
  u_graph::CFG g = g11();
  if (DEBUG) {
    std::cout << "Input Graph" << "\n\n";
    g.print_dot();
  }

  spanning_tree::Tree t = g.compute_spanning_tree();
  if (DEBUG) {
    std::cout << "\n\n" << "Spanning tree" << "\n\n";
    t.print_dot();
  }

  vst::cycle_equiv(t);
  if (DEBUG) {
    std::cout << "\n\n" << "Updated Spanning tree" << "\n\n";
    t.print_dot();
  }

  vst::VST tv = vst::VST(t);
  if (DEBUG) {
    std::cout << "\n\n" << "VST" << "\n\n";
    tv.print_dot();
  }

  // pst::PST p = pst::PST(t);
  // if (DEBUG) {
  //   std::cout << "\n\n" << "PST" << "\n\n";
  //   std::cout << "PST size: " << p.size() << "\n";
  //   p.print_dot();
  // }

  tree::Tree p = pst::compute_pst(t);
  if (DEBUG) {
    std::cout << "\n\n" << "PST" << "\n\n";
    p.print_dot();
  }


  tree::Tree pv = vst::compute_pvst(t);
  if (DEBUG) {
    std::cout << "\n\n" << "PVST" << "\n\n";
    pv.print_dot();
  }
  */

  return  0;
}
