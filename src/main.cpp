#include <iostream>

#include "./graph/graph.hpp"
#include "./graph/spanning_tree.hpp"
#include "./vst/vst.hpp"

int main() {
  //u_graph::CFG k = u_graph::CFG::CFG();
  u_graph::CFG g;

  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(2, 3);
  g.add_edge(3, 1);

  g.set_start_node(0);
  g.set_stop_node(3);

  //g.print_dot();

  spanning_tree::Tree t = g.compute_spanning_tree();

  // t.print_dot();

  vst::cycle_equiv(t);

  t.print_dot();


  vst::VST tv = vst::VST(t);

  tv.print_dot();

  return  0;
  g = u_graph::CFG();

  //g.print_dot();

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
  g.set_stop_node(9);
  g.print_dot();

  t = g.compute_spanning_tree();

  t.print_dot();
}
