#ifndef BIDIRECTED_HPP
#define BIDIRECTED_HPP


#include <cstddef>
#include <string>
#include <vector>
#include <unordered_set>
#include <set>

#include <handlegraph/handle_graph.hpp>
//#include <handlegraph/mutable_handle_graph.hpp>

namespace hg = handlegraph;



namespace bidirected {

/**
 * Path (or color)
 * ---------------
 */
struct path_t {
  std::string name; // name as pertains the GFA file
  std::size_t id; // numerical id to be associated with handle
  bool is_circular; // is the path circular?
};


/**
   l (left) 5' or +
   r (left) 3' or -
 */
enum class VertexEnd {
	l,
	r
};



struct PathInfo {
  std::size_t path_id;
  std::size_t step_index;

  PathInfo(): path_id(0), step_index(0) {}
  PathInfo(std::size_t path_id, std::size_t step_index): path_id(path_id), step_index(step_index) {}
};

/**
   Each edge connects a pair of vertices.
   An edge is + incident or - incident to each vertex which is a VertexEnd
   the pair is vertex id and vertex incident side
   all edges in this graph are gray and have no labels
 */
class Edge {
  std::size_t v1_idx;
  VertexEnd v1_end;
  std::size_t v2_idx;
  VertexEnd v2_end;

  //std::pair<std::size_t, VertexEnd> v1;
  //std::pair<std::size_t, VertexEnd> v2;

public:
  // constructor
  Edge();
  Edge(std::size_t v1, VertexEnd v1_end, std::size_t v2, VertexEnd v2_end);

  // getters
  std::size_t get_v1_idx() const;
  VertexEnd get_v1_end() const;
  std::size_t get_v2_idx() const;
  VertexEnd get_v2_end() const;

  // operators
  // << operator
  friend std::ostream& operator<<(std::ostream& os, const Edge& e);
};

/**
   a vertex is invalid if it lacks either a label or a handle
 */
class Vertex {
  std::string label; // or sequence
  //std::unordered_set<std::size_t> edges; // indexes to the edge vector in Graph

  // indexes to the edge vector in Graph
  std::set<std::size_t> edges_l;
  std::set<std::size_t> edges_r;

  // paths (also colors)
  // the first element is the path id and the second is the step index or the coordinate of the sequence in that linear haplotype
  // TODO: change to unordered_set
  std::vector<PathInfo> paths;

  // from libHandleGraph
  std::string handle;
  bool is_reversed_;

public:
  // ------------
  // constructors
  // ------------
  Vertex();
  Vertex(const std::string& label);
  Vertex(const std::string& label, const handlegraph::nid_t& id);
  //Vertex(std::string label, std::unordered_set<std::size_t> edge_index);

  /*
	setters and getters
	-------------------
  */

  // -------
  // getters
  // -------

  bool is_reversed() const;
  const std::string& get_label() const;
  const std::string& get_handle() const;
  const std::set<std::size_t>& get_edges_l() const;
  const std::set<std::size_t>& get_edges_r() const;

  // -------
  // setters
  // -------

  // returns the new value
  bool toggle_reversed();
  //void set_handle(const std::string& handle);
  //void set_label(const std::string& label);
  void add_edge(std::size_t edge_index, VertexEnd vertex_end);
  //void remove_edge(std::size_t edge_index);
  //std::unordered_set<std::size_t> get_edges() const;

  // It is up to the user to make sure that the path_id is not already in the "set"
  void add_path(std::size_t path_id, std::size_t step_index);

};

/**
   A variation graph as a bidirected graph
 */
class VariationGraph {
  std::vector<Vertex> vertices;
  std::vector<Edge> edges;

  std::vector<path_t> paths;

  std::unordered_set<std::size_t> start_nodes;
  std::unordered_set<std::size_t> end_nodes;

  // for libHandleGraph
  // min and max vertex ids
  std::size_t min_id;
  std::size_t max_id;

public:
  // -----------
  // constructor
  // -----------
  VariationGraph();
  VariationGraph(std::size_t vertex_count, std::size_t edge_count, std::size_t path_count);


  // -------
  // getters
  // -------
  std::size_t size() const; // the number of vertices in the graph valid or not
  const Vertex& get_vertex(std::size_t index) const;
  Vertex& get_vertex_mut(std::size_t index);
  const Edge& get_edge(std::size_t index) const;


  // -------
  // setters
  // --------
  void append_vertex();   // adds an invalid vertex to the graph
  void add_vertex(const Vertex& vertex);
  void add_edge(std::size_t v1, VertexEnd v1_end, std::size_t v2, VertexEnd v2_end);

  void add_start_node(std::size_t node_id);
  void add_stop_node(std::size_t node_id);

  void set_min_id(std::size_t min_id);
  void set_max_id(std::size_t max_id);


  // ----
  // misc
  // ----
  void dbg_print();


  /*
	Library implementations
	=======================

	- HandleGraph
	- MutableHandleGraph
	- MutablePathHandleGraph
  */

  /*
	 HandleGraph
	 -----------

	Partially implements HandleGraph

	Implemented:
	 -
   */


  /// Method to check if a node exists by ID
  bool has_node(handlegraph::nid_t node_id) const;

  /// Look up the handle for the node with the given ID in the given orientation
  // returns the node_id as a handle
  // TODO: not applicable, make applidable
  handlegraph::handle_t get_handle(const handlegraph::nid_t& node_id,
								   bool is_reverse = false) const;

  /// Get the ID from a handle
  // same as a size_t cast
  handlegraph::nid_t get_id(const handlegraph::handle_t& handle) const;

  /// Get the orientation of a handle
  bool get_is_reverse(const handlegraph::handle_t& handle) const;

  /// Invert the orientation of a handle (potentially without getting its ID)
  /*
	since vertex doesn't have orientation yet this will change the incidence of the edges
  edges incident with the l side will be incident with the r side and vice versa
   */
  handlegraph::handle_t flip(const handlegraph::handle_t& handle);

  /// Get the length of a node (label in bp)
  size_t get_length(const handlegraph::handle_t& handle) const;

  /// Get the sequence of a node, presented in the handle's local forward
  /// orientation.
  std::string get_sequence(const handlegraph::handle_t& handle) const;

  /// Return the number of nodes in the graph
  std::size_t get_node_count() const;

  /// Return the smallest ID in the graph, or some smaller number if the
  /// smallest ID is unavailable. Return value is unspecified if the graph is empty.
  handlegraph::nid_t min_node_id() const;

  /// Return the largest ID in the graph, or some larger number if the
  /// largest ID is unavailable. Return value is unspecified if the graph is empty.
  handlegraph::nid_t max_node_id() const;

  // TODO: implement
  bool follow_edges_impl(const handlegraph::handle_t& handle,
						 bool go_left,
						 const std::function<bool(const handlegraph::handle_t&)>& iteratee) const;


  // TODO: implement
  bool for_each_handle_impl(const std::function<bool(const handlegraph::handle_t&)>& iteratee,
							bool parallel = false) const;


	/*
	  MutableHandleGraph
	  -------------------

	partially implements MutableHandleGraph

	Implemented:
	 - handlegraph::handle_t create_handle(const std::string& sequence);
	 - handlegraph::handle_t create_handle(const std::string& sequence, const handlegraph::nid_t& id);
   */


  // adds a sequence to the end of the graph
  handlegraph::handle_t create_handle(const std::string& sequence);

  handlegraph::handle_t create_handle(const std::string& sequence,
									  const handlegraph::nid_t& id);

  // this seems to assume a biedged representation
  void create_edge(const handlegraph::handle_t& left, const handlegraph::handle_t& right);

  /*
	MutablePathHandleGraph
	----------------------

	partially implements MutablePathHandleGraph

	Implemented:
	 - create_path_handle
   */

   /**
	 * Create a path with the given name. The caller must ensure that no path
	 * with the given name exists already, or the behavior is undefined.
	 * Returns a handle to the created empty path. Handles to other paths must
	 * remain valid.
	 */
  hg::path_handle_t create_path_handle(const std::string& name,
									   bool is_circular = false);


  /**
   * Renames a path. Existing path_handle_t's may become invalidated..
   */
  hg::path_handle_t rename_path(const hg::path_handle_t& path_handle,
								const std::string& new_name);

	/**
	 * Append a visit to a node to the given path. Returns a handle to the new
	 * final step on the path which is appended. If the path is cirular, the new
	 * step is placed between the steps considered "last" and "first" by the
	 * method path_begin. Handles to prior steps on the path, and to other paths,
	 * must remain valid.
	 */
  hg::step_handle_t append_step(const hg::path_handle_t& path,
								const hg::handle_t& to_append);

  /**
	 * Make a path circular or non-circular. If the path is becoming circular, the
	 * last step is joined to the first step. If the path is becoming linear, the
	 * step considered "last" is unjoined from the step considered "first" according
	 * to the method path_begin.
	 */
  void set_circularity(const hg::path_handle_t& path, bool circular);

};


}; // namespace bidirected
#endif
