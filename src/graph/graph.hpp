#ifndef GRAPH_HPP
#define GRAPH_HPP

#include <cstddef>
#include <vector>
#include <set>

#include "./spanning_tree.hpp"
#include "./digraph.hpp"

/*
 * TODO:
 * - rename namespace to flow graph
 * - use smart pointers
 */

// undirected graph
namespace u_graph {

/**
 * An edge in an undirected graph where (l , r) == (r, l)
 */
class Edge {
  std::size_t l;
  std::size_t r;
public:
  Edge(std::size_t l, std::size_t r):  l(std::min(l,r)), r(std::max(l,r)){}

// spaceship operator
friend constexpr auto operator<=>(Edge, Edge) = default;
};

class Vertex {
  std::set<std::size_t> adj_vertices;
public:
  Vertex();

  void add_adjacent_vertex(std::size_t vertex);
  void del_adjacent_vertex(std::size_t vertex);

  std::set<std::size_t>const& get_adjacent_vertices() const;
};

/**
 * this undirected graph is actually a flow graph
 * undirected
 * edge from end to start will always be a backedge
 */
class CFG {
  std::vector<Vertex> adj_list;
  // start node is always zero
  // stop node is always N+1
  // std::set<std::size_t> start_nodes;
  // std::set<std::size_t> end_nodes;
  std::size_t start_node_id{}; // static data member

  /*
    private methods
   */
  std::size_t stop_node_internal_idx();
  std::size_t size_internal();

  Vertex& start_node_internal();
  Vertex& stop_node_internal();
  Vertex const& get_vertex_internal(std::size_t vertex) const;
  std::set<std::size_t> const& get_adjacent_vertices_internal(std::size_t vertex) const;

public:
  // CFG();
  // TODO: from gfa
  // CFG(std::size_t initial_len=2); // from di graph or from gfa
  CFG(std::size_t initial_len=2);
  CFG(digraph::DiGraph const& di_graph);


  // setters
  void add_edge(std::size_t n1, std::size_t n2);

  Vertex& get_vertex_mut(std::size_t vertex);
  Vertex const& get_vertex(std::size_t vertex) const;
  std::set<std::size_t> const& get_adjacent_vertices(std::size_t vertex) const;

  // mark/set a vertex as a start node
  void set_start_node(std::size_t vertex);
  void set_stop_node(std::size_t vertex);
  std::size_t size();
  spanning_tree::Tree compute_spanning_tree();

  // split any node with > 1 incoming and >1 outgoing edges
  // into two nodes connected by a single (grey) edge
  void make_bi_edged();

  void print_dot();
};

}; // namespace u_graph
#endif
