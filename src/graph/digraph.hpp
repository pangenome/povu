#ifndef DIGRAPH_HPP
#define DIGRAPH_HPP

#include <cstddef>
#include <iostream>
#include <unordered_set>
#include <set>
#include <vector>

namespace digraph {

enum colour { grey, black };
/*
 * Edge
 * ----
 */
class Edge {
  std::size_t frm; // from
  std::size_t t; // to
  colour c;

public:
  Edge();
  Edge(std::size_t frm, std::size_t to, colour c=colour::black);

  std::size_t to() const;
  std::size_t from() const;
  colour get_colour() const;

  void set_from(std::size_t f);
  void set_to(std::size_t t);
};

bool operator<(const Edge& lhs, const Edge& rhs);

/*
 * Vertex
 * ------
 */
// TODO: use unordered set?
class Vertex {
  std::set<Edge> o;
  std::set<Edge> i;

public:
  Vertex();

  std::set<Edge> const& out() const;
  std::set<Edge> const& in() const;


  std::set<Edge>* out_mut ();
  std::set<Edge>* in_mut();

  
  void add_out(std::size_t self_idx, std::size_t to_idx);
  void add_in(std::size_t from_idx, std::size_t self_idx);

  bool is_leaf() const;
};

/*
 * Vertex
 * ------
 */
class DiGraph {
  std::vector<Vertex> adj; // adjacency list to store edges
  std::set<std::size_t> start_nodes;
  std::set<std::size_t> end_nodes;

public:
  DiGraph();
  DiGraph(std::size_t size);
  DiGraph(std::set<std::size_t>&& start_nodes, std::set<std::size_t>&& stop_nodes);

  void add_start_node(std::size_t idx);
  void add_stop_node(std::size_t idx);

  // tree::Tree spanning_tree();
  std::size_t size() const;

  Vertex const& get_vertex(std::size_t idx) const;
  Vertex& get_vertex_mut(std::size_t idx);


  std::set<std::size_t> const &starts() const;
  std::set<std::size_t> const &stops() const;

  void add_edge(std::size_t from, std::size_t to);

  // convert the digraph into a biedged graph
  // nodes with >1 incoming edges and >1 outgoing edges
  // are split 2 nodes connected by a grey directed edge
  void biedge();

  void print_dot();

};
} // namespace digraph
#endif
