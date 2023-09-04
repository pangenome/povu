#include <iostream>

#include "./graph/digraph.hpp"
#include "./graph/u_graph.hpp"
#include "./graph/spanning_tree.hpp"
#include "./vst/vst.hpp"
#include "./vst/pst.hpp"
#include "./pvst/pvst.hpp"
#include "graph/tree.hpp"


#include <handlegraph/handle_graph.hpp>
#include <handlegraph/mutable_handle_graph.hpp>
#include <handlegraph/mutable_path_deletable_handle_graph.hpp>

u_graph::FlowGraph g1() {
  u_graph::FlowGraph g;

  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(2, 3);
  g.add_edge(3, 1);

  g.set_start_node(0);
  g.set_stop_node(3);

  return g;
};

// big
u_graph::FlowGraph g2() {
  u_graph::FlowGraph g;
  g = u_graph::FlowGraph();

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
u_graph::FlowGraph g3() {
  u_graph::FlowGraph g;
  g = u_graph::FlowGraph();

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
u_graph::FlowGraph g4() {
  u_graph::FlowGraph g;
  g = u_graph::FlowGraph();

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
u_graph::FlowGraph g5() {
  u_graph::FlowGraph g;
  g = u_graph::FlowGraph();

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
u_graph::FlowGraph g6() {
  u_graph::FlowGraph g;

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

digraph::DiGraph g13() {
  digraph::DiGraph g;

  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(2, 3);
  g.add_edge(3, 1);
  g.add_edge(3, 4);

  g.add_start_node(0);
  g.add_stop_node(4);

  return g;
};

digraph::DiGraph g14() {
  digraph::DiGraph g;

  g.add_edge(0, 1);

  g.add_edge(1, 2);
  g.add_edge(1, 3);

  g.add_edge(2, 4);
  g.add_edge(3, 4);

  g.add_edge(4, 5);

  g.add_start_node(0);
  g.add_stop_node(5);

  return g;
};

digraph::DiGraph g15() {
  digraph::DiGraph g;

  g.add_edge(0, 1);

  g.add_edge(1, 2);
  g.add_edge(1, 3);

  g.add_edge(2, 4);
  g.add_edge(3, 4);

  g.add_edge(4, 5);
    g.add_edge(1, 5);
    g.add_edge(5, 6);


  g.add_start_node(0);
  g.add_stop_node(6);

  return g;
};

// figure 1 Johnson et al
digraph::DiGraph g16() {
  digraph::DiGraph g;

  g.add_edge(0, 1);

  //g.add_edge(0, 1);
  g.add_edge(1, 3);
  g.add_edge(3, 5);
  g.add_edge(5, 8);
  g.add_edge(5, 9);
  g.add_edge(8, 10);
  g.add_edge(9, 10);
  g.add_edge(10, 11);
  g.add_edge(11, 1);
  g.add_edge(11, 15);
  g.add_edge(15, 16);

  g.add_edge(0, 2);
  g.add_edge(2, 4);
  g.add_edge(4, 6);
  g.add_edge(4, 7);
  g.add_edge(6, 12);
  g.add_edge(7, 12);
  g.add_edge(12, 13);
  g.add_edge(13, 7);
  g.add_edge(13, 14);
  g.add_edge(14, 15);
  g.add_edge(15, 16);

  g.add_start_node(0);
  g.add_stop_node(16);

  return g;
};

// diamond
digraph::DiGraph g17() {
  digraph::DiGraph g;
  //g = u_graph::FlowGraph();

  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(1, 3);
  g.add_edge(2, 4);
  g.add_edge(3, 4);
  g.add_edge(4, 5);

  g.add_start_node(0);
  g.add_stop_node(5);

  return g;
}

// with overlapping
digraph::DiGraph g18() {
  digraph::DiGraph g;

  g.add_edge(0, 1);
  g.add_edge(1, 3);
  g.add_edge(1, 2);
  g.add_edge(2, 4);

  g.add_edge(3, 4);
  g.add_edge(5, 3);

  g.add_edge(4, 5);
  g.add_edge(5, 6);

  g.add_start_node(0);
  g.add_stop_node(6);

  return g;
};

// with capping
digraph::DiGraph g19() {
  digraph::DiGraph g;

  // TODO: flip 2 3 for capping
  g.add_edge(0, 1);
  g.add_edge(1, 3);
  g.add_edge(1, 2);
  g.add_edge(2, 4);

  g.add_edge(3, 4);
  g.add_edge(5, 2);

  g.add_edge(4, 5);
  //g.add_edge(5, 2);

  g.add_edge(5, 6);

  g.add_start_node(0);
  g.add_stop_node(6);

  return g;
};

// figure 1 Johnson et al
digraph::DiGraph g20() {
  digraph::DiGraph g;

  g.add_edge(0, 1);

  //g.add_edge(0, 1);
  g.add_edge(1, 3);
  g.add_edge(3, 5);
  g.add_edge(5, 8);
  g.add_edge(5, 9);
  g.add_edge(8, 10);
  g.add_edge(9, 10);
  g.add_edge(10, 11);
  g.add_edge(11, 1);
  g.add_edge(11, 15);
  g.add_edge(15, 16);

  g.add_edge(0, 2);
  g.add_edge(2, 4);
  g.add_edge(4, 6);
  g.add_edge(4, 7);
  g.add_edge(6, 12);
  g.add_edge(7, 12);
  g.add_edge(12, 13);
  g.add_edge(13, 7);
  g.add_edge(13, 14);
  g.add_edge(14, 15);
  g.add_edge(15, 16);

  g.add_start_node(0);
  g.add_stop_node(16);

  return g;
};


bool DEBUG = false;

int main() {
  
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

  vst::cycle_equiv(t);
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
