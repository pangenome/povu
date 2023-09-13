#ifndef DIGRAPH_HPP
#define DIGRAPH_HPP

#include <cstddef>
#include <iostream>
#include <string>
#include <unordered_set>
#include <set>
#include <vector>
#include <functional>

#include "../core/core.hpp"

#include <handlegraph/handle_graph.hpp>
//#include <handlegraph/mutable_handle_graph.hpp>

namespace hg = handlegraph;

namespace digraph {

/*
 * Path
 * ----
 */
struct path_t {
  std::string name; // name as pertains the GFA file
  std::size_t id; // numerical id to be associated with handle
};

bool operator<(const path_t& lhs, const path_t& rhs);
  
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
  std::string seq;
  // todo: store a size_t ?
  std::string handle; // name, id or handle lack of a handle means the vertex is invalid/unset

  // path ids
  // a path should be a more complex object (see odgi)
  // name of the path the vertex is in
  std::set<std::size_t> paths;

  // TODO: add orientation

public:
  Vertex();
  Vertex(const std::string& sequence);
  Vertex(const std::string& sequence, const handlegraph::nid_t& id);

  std::set<Edge> const& out() const;
  std::set<Edge> const& in() const;

  std::string const& get_seq() const;
  std::string const& get_handle() const;


  std::set<Edge>* out_mut ();
  std::set<Edge>* in_mut();

  void add_out(std::size_t self_idx, std::size_t to_idx, core::color c);
  void add_in(std::size_t from_idx, std::size_t self_idx, core::color c);

  //void set_seq(std::string& s);
  //void set_handle(std::string&& h);

  // returns a positive value when the vertex is valid and the handle is set and
  // path was empty
  int set_path(std::size_t p_id);


  bool is_leaf() const;
};

/*
 * Vertex
 * ------
 */

// TODO: handle the case of disconnected components
class DiGraph : public handlegraph::HandleGraph {
  std::vector<Vertex> adj; // adjacency list to store edges

  // nodes which have no incoming edges and no outgoing edges do are not yet handled/stored
  std::set<std::size_t> start_nodes; // nodes with no incoming edges
  std::set<std::size_t> end_nodes; // nodes with no outgoing edges

  std::vector<path_t> paths; // names of the paths in the graph, does not enumerate the nodes in the paths

public:

  // constructor(s)
  // --------------
  DiGraph();
  DiGraph(std::size_t size);
  DiGraph(std::size_t size, std::size_t path_count);
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

  void compute_start_nodes();
  void compute_stop_nodes();

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

  /*
	Library implementations
  */

  // HandleGraph
  // -----------

  bool has_node(handlegraph::nid_t node_id) const;

  // returns the node_id as a handle
  // TODO: not applicable, make applidable
  handlegraph::handle_t get_handle(const handlegraph::nid_t& node_id,
								   bool is_reverse = false) const;

  // same as a size_t cast
  handlegraph::nid_t get_id(const handlegraph::handle_t& handle) const;

  // FIXME: always returns false for now
  bool get_is_reverse(const handlegraph::handle_t& handle) const;

  // since vertex doesn't have orientation yet this is useless and just returns the handle
  handlegraph::handle_t flip(const handlegraph::handle_t& handle) const;

  /// Get the length of a node
  // TODO: not applicable, make applidable
  size_t get_length(const handlegraph::handle_t& handle) const;

  std::string get_sequence(const handlegraph::handle_t& handle) const;

  std::size_t get_node_count() const;

  handlegraph::nid_t min_node_id() const;

  handlegraph::nid_t max_node_id() const;

 // TODO: implement
 bool follow_edges_impl(const handlegraph::handle_t& handle,
						bool go_left,
						const std::function<bool(const handlegraph::handle_t&)>& iteratee) const;


  // TODO: implement
  bool for_each_handle_impl(const std::function<bool(const handlegraph::handle_t&)>& iteratee,
							bool parallel = false) const;


  // MutableHandleGraph
  // ------------------
  // TODO: not yet fully implements MutableHandleGraph

  // adds a sequence to the end of the graph
  handlegraph::handle_t create_handle(const std::string& sequence);

  handlegraph::handle_t create_handle(const std::string& sequence,
									  const handlegraph::nid_t& id);

  // the second arg is expected to be a string castable into size_t
  //handlegraph::handle_t create_handle(const std::string& sequence, const std::string const& name);

  void create_edge(const handlegraph::handle_t& left,
				   const handlegraph::handle_t& right);

  // MutablePathHandleGraph
  // ----------------------

  // TODO: not yet fully implements MutablePathHandleGraph
  // TODO: is circular is not implemented


  hg::path_handle_t create_path_handle(const std::string& name,
									   bool is_circular = false);


  hg::path_handle_t rename_path(const hg::path_handle_t& path_handle,
								const std::string& new_name);

};
} // namespace digraph
#endif
