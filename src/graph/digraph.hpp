#ifndef DIGRAPH_HPP
#define DIGRAPH_HPP

#include <cstddef>
#include <iostream>
#include <unordered_set>
#include <set>
#include <vector>

#include "../core/core.hpp"

#include <handlegraph/mutable_handle_graph.hpp>

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

// TODO: handle the case of disconnected components
class DiGraph {
  std::vector<Vertex> adj; // adjacency list to store edges
  std::set<std::size_t> start_nodes;
  std::set<std::size_t> end_nodes;

public:

  //
  // ------------

  // add a sequence to the end of the graph
  // increment the size of the graph by 1
  handlegraph::handle_t create_handle(const std::string& sequence);

  handlegraph::handle_t create_handle(const std::string& sequence,
									  const handlegraph::nid_t& id);


  /// Create an edge connecting the given handles in the given order and orientations.
  /// Ignores existing edges.
  void create_edge(const handlegraph::handle_t& left,
				   const handlegraph::handle_t& right);


  handlegraph::handle_t apply_orientation(const handlegraph::handle_t& handle);

  std::vector<handlegraph::handle_t>
  divide_handle(const handlegraph::handle_t& handle, const std::vector<std::size_t>& offsets);

  void optimize(bool allow_id_reassignment = true);


  bool apply_ordering(const std::vector<handlegraph::handle_t>& order, bool compact_ids = false);

  void set_id_increment(const handlegraph::nid_t& min_id);

  void increment_node_ids(handlegraph::nid_t increment);
  
  void increment_node_ids(long increment);
    
  void reassign_node_ids(const std::function<handlegraph::nid_t(const handlegraph::nid_t&)>& get_new_id);
  
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

  // if the from node or to node does not exist, add edges until they
  // do exist
  
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
