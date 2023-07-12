#ifndef DIGRAPH_HPP
#define DIGRAPH_HPP

#include <cstddef>
#include <iostream>
#include <unordered_set>
#include <set>
#include <vector>

#include "../core/core.hpp"

namespace digraph {

/*
 * Edge
 * ----
 */
class Edge {
  std::size_t frm; // from
  std::size_t t; // to
  core::color c;

public:
  Edge();
  Edge(std::size_t frm, std::size_t to, core::color c=core::color::black);

  std::size_t to() const;
  std::size_t from() const;

  core::color get_color() const;
  bool is_black() const;

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

  void add_out(std::size_t self_idx, std::size_t to_idx, core::color c);
  void add_in(std::size_t from_idx, std::size_t self_idx, core::color c);

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
  // constructor(s)
  // --------------
  DiGraph();
  DiGraph(std::size_t size);
  DiGraph(std::set<std::size_t>&& start_nodes, std::set<std::size_t>&& stop_nodes);


  // getters
  // -------
  std::set<std::size_t> const &starts() const;
  std::set<std::size_t> const &stops() const;

  // tree::Tree spanning_tree();
  std::size_t size() const;

  Vertex const& get_vertex(std::size_t idx) const;

  // setters
  // -------
  void add_start_node(std::size_t idx);
  void add_stop_node(std::size_t idx);

  Vertex& get_vertex_mut(std::size_t idx);

  /*
   * add an edge to the graph
   * if the edge already exists, do nothing
   * if the edge is a loop, still add it
   * if the edge is a loop, and the color is grey, do nothing FIXME??
   */
  void add_edge(std::size_t from, std::size_t to, core::color c=core::color::black);

  /*
   * convert the digraph into a biedged graph
   * nodes with >1 incoming edges or >1 outgoing edges
   * are split 2 nodes connected by a grey directed edge
   */
  void biedge();

  void print_dot();
};
} // namespace digraph
#endif
