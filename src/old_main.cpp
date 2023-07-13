#include <iostream>
#include "./graph/graph.hpp"
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
  g.print_dot();

  spanning_tree::Tree t = g.compute_spanning_tree();

  t.print_dot();


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

int old_main() {

  graph::Graph g = graph::Graph();
  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(2, 3);
  g.add_edge(3, 1);

  g.add_start_node(0);
  g.add_stop_node(3);
  g.print_dot();

  graph::FlowGraph f = graph::FlowGraph(g);
  f.print_dot();

  tree::Tree t = f.spanning_tree();
  t.print_dot();

  // --------------
  g = graph::Graph();
  g.add_edge(0, 1);
  g.add_edge(1, 2);
  g.add_edge(2, 1);
  g.add_edge(2, 3);
  g.add_edge(3, 0);
  g.add_edge(3, 4);
  g.add_edge(4, 4);
  g.print_dot();

  g.add_start_node(0);
  g.add_stop_node(4);

  f = graph::FlowGraph(g);
  f.print_dot();

  t = f.spanning_tree();
  t.print_dot();

  // --------------
  g = graph::Graph();
  g.add_edge(0, 1);
  g.add_edge(0, 2);
  g.add_edge(1, 3);
  g.add_edge(2, 3);
  g.add_edge(3, 4);
  g.add_edge(0, 4);

  g.add_start_node(0);
  g.add_stop_node(4);

  g.print_dot();

  f = graph::FlowGraph(g);
  f.print_dot();

  t = f.spanning_tree();
  t.print_dot();

  return 0;
}
