#ifndef BIDIRECTED_HPP
#define BIDIRECTED_HPP

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <sys/types.h>
#include <vector>

#include <handlegraph/handle_graph.hpp>
#include "../core/core.hpp"
#include "../core/constants.hpp"
#include "../common/common.hpp"

namespace hg = handlegraph;

namespace bidirected {

using namespace graph_types;


// TODO: replace with struct or class to allow methods like complement
enum class orientation_t {
  forward,
  reverse
};
std::ostream& operator<<(std::ostream& os, const orientation_t& o);

struct id_n_orientation_t {
  std::size_t v_idx;
  orientation_t orientation;
};
std::ostream& operator<<(std::ostream& os, const id_n_orientation_t& x);


struct PathInfo {
  std::size_t path_id;
  std::size_t step_index;

  // --------------
  // constructor(s)
  // --------------
  PathInfo(): path_id(0), step_index(0) {}
  PathInfo(std::size_t path_id, std::size_t step_index): path_id(path_id), step_index(step_index) {}
};

/**
 * Each edge connects a pair of vertices.
 * An edge is + incident or - incident to each vertex which is a VertexEnd
 * the pair is vertex id and vertex incident side
 * all edges in this graph are gray and have no labels
 */
class Edge {
  std::size_t v1_idx;
  VertexEnd v1_end;
  std::size_t v2_idx;
  VertexEnd v2_end;

  // the equivalence class of the edge
  std::size_t eq_class {core::constants::UNDEFINED_SIZE_T};

public:
  // --------------
  // constructor(s)
  // --------------
  Edge();
  Edge(std::size_t v1, VertexEnd v1_end, std::size_t v2, VertexEnd v2_end);

  // ---------
  // getter(s)
  // ---------
  std::size_t get_v1_idx() const;
  VertexEnd get_v1_end() const;
  std::size_t get_v2_idx() const;
  VertexEnd get_v2_end() const;
  // get the other vertex and side connected to this edge
  // does not check to confirm that the vertex is connected to this edge
  // TODO:
  //   - nice to have:
  //     * check that the vertex at vertex index is connected to this edge
  side_n_id_t get_other_vertex(std::size_t vertex_index) const;
  std::size_t get_eq_class() const;

  // ---------
  // setter(s)
  // ---------
  void set_v1_idx(std::size_t v1_idx);
  // void set_v1_end(VertexEnd v1_end);
  void set_v2_idx(std::size_t v2_idx);
  // void set_v2_end(VertexEnd v2_end);
  void set_eq_class(std::size_t eq_class);

  // -----------
  // operator(s)
  // -----------
  friend std::ostream& operator<<(std::ostream& os, const Edge& e);
};


/**
 * a vertex is invalid if it lacks either a label or a handle
 */
class Vertex {
  std::string label; // or sequence

  // indexes to the edge vector in Graph
  std::set<std::size_t> edges_l;
  std::set<std::size_t> edges_r;

  // paths (also colors)
  // the first element is the path id and the second is the step index or
  // the coordinate of the sequence in that linear haplotype
  // TODO: change to unordered_set
  std::vector<PathInfo> paths;

  // from libHandleGraph
  std::string handle;
  bool is_reversed_;

  std::string name_; // sequence name in the GFA file

  // the equivalence class of the vertex
  std::size_t eq_class {core::constants::UNDEFINED_SIZE_T};

public:
  // --------------
  // constructor(s)
  // --------------
  Vertex();
  Vertex(const std::string& label);
  Vertex(const std::string& label, const handlegraph::nid_t& id);


  // ---------
  // getter(s)
  // ---------
  bool is_reversed() const;
  const std::string& get_label() const;
  // rc is reverse complement
  std::string get_rc_label() const;
  const std::string& get_handle() const;
  const std::string& get_name() const;
  const std::set<std::size_t>& get_edges_l() const;
  const std::set<std::size_t>& get_edges_r() const;
  const std::vector<PathInfo>& get_paths() const;
  std::size_t get_eq_class() const;


  // ---------
  // setter(s)
  // ---------
  // returns the new value
  bool toggle_reversed();
  void add_edge(std::size_t edge_index, VertexEnd vertex_end);
  void clear_edges();
  // void set_label(const std::string& label);
  void set_handle(const std::string& handle);
  void set_handle(id_t id);
  void set_name(const std::string& name);

  // It is up to the user to make sure that the path_id is not already in the "set"
  void add_path(std::size_t path_id, std::size_t step_index);
  void set_eq_class(std::size_t eq_class);

};

struct component;



/**
 * A variation graph as a bidirected graph
 */
class VariationGraph {
  std::vector<Vertex> vertices;
  std::vector<Edge> edges;

  // TODO: make this a map for easier work with components
  // we want unique path IDs for the entire graph
  std::map<id_t, path_t> paths;

  // TODO: associate with paths above, maybe make it a map as well or merge them into one
  std::vector<std::vector<id_n_orientation_t>> raw_paths;

  // vertices with incident edges on only one side
  std::set<side_n_id_t> tips_;

  // graph start nodes are vertices with edges adjacent to only one side, the right side or the left side
  // the side is the one without any incident edges
  //std::set<std::size_t> graph_start_nodes_;
  //std::set<std::size_t> graph_end_nodes_;

  // we store the side which would visit a black edge
  // haplotype start nodes are vertices which start paths according to the P lines in a GFA file
  // start side is the side from which the path starts
  // e.g. 5+ will have a start side of 5 left and 5- will have a start side of 5 right
  std::set<side_n_id_t> haplotype_start_nodes_;
  std::set<side_n_id_t> haplotype_end_nodes_;

  // for libHandleGraph
  // min and max vertex ids
  std::size_t min_id;
  std::size_t max_id;

public:
  // --------------
  // constructor(s)
  // --------------
  VariationGraph();
  VariationGraph(std::size_t vertex_count, std::size_t edge_count, std::size_t path_count);


  // ---------
  // getter(s)
  // ---------
  std::size_t size() const; // the number of vertices in the graph valid or not
  const Vertex& get_vertex(std::size_t index) const;

  const Vertex& get_vertex_by_name(std::string n) const;
  std::size_t get_vertex_idx_by_name(std::string n) const;


  Vertex& get_vertex_mut(std::size_t index);
  const Edge& get_edge(std::size_t index) const;
  Edge& get_edge_mut(std::size_t index);
  const std::vector<Edge>& get_all_edges() const;

  /**
   * @brief Get all paths between two nodes in the bidirected variation graph
   *
   * @param start_id the start node id
   * @param stop_id the stop node id
   * @param compact if true, remove redundant strand information in the path
   * @return std::vector<subpaths_t> a vector of paths
   */
  std::vector<std::vector<side_n_id_t>> get_paths(id_t start_id, id_t stop_id, bool compact=true) const;

  // get adjacent vertex indexes to a vertex in a given direction
  std::vector<side_n_id_t> get_adj_vertices(std::size_t vertex_index, VertexEnd vertex_end) const;
  //std::unordered_set<std::size_t> get_start_nodes() const;
  //std::unordered_set<std::size_t> get_end_nodes() const;
  std::vector<path_t> get_paths() const; // TODO: rename to get_haplotypes // return a reference to the map
  const path_t& get_path(std::size_t path_id) const;
  std::size_t get_path_count() const;

  std::set<side_n_id_t> const& tips() const;

  // tips that are not hap starts or hap ends
  // by definition they are not part of any path as well
  std::set<side_n_id_t> get_orphan_tips() const;

  // graph start nodes are tips that are haplotype starts or ends
  // if strict is true, only return tips which are haplotype start or end nodes
  // i.e. non orphan tips
  // TODO: implement strict
  std::set<side_n_id_t> graph_start_nodes(bool strict=false) const;

  // if strict is false it allows orphan tips
  std::set<side_n_id_t> graph_end_nodes(bool strict=false) const;

  // TODO: remove the find part? replace with get?
  std::set<side_n_id_t> const& find_haplotype_start_nodes() const;
  std::set<side_n_id_t> const& find_haplotype_end_nodes() const;

  // TODO: maybe nice to have?
  //std::set<std::size_t> get_edges(std::size_t vertex_index, VertexEnd vertex_end) const;

  /**
   * @brief Get the minimum vertex id
   */
  bool validate_haplotype_paths();

  // ---------------------
  // setter(s) & modifiers
  // ---------------------
  void append_vertex();   // adds an invalid vertex to the graph
  /**
    * @brief Add a vertex to the graph
    *
    * @param vertex the vertex to add
    * @return std::size_t the index of the vertex in the graph
   */
  std::size_t add_vertex(const Vertex& vertex);


  void add_edge(const Edge& edge); // handles the case where one or both of the vertices are invalid
  std::size_t add_edge(std::size_t v1, VertexEnd v1_end, std::size_t v2, VertexEnd v2_end);

  void add_path(const path_t &path);
  void set_raw_paths(std::vector<std::vector<id_n_orientation_t>> &raw_paths);

  void add_tip(std::size_t node_id, VertexEnd end);

  //bool add_graph_start_node(std::size_t node_id);
  //bool add_graph_end_node(std::size_t node_id);

  void add_haplotype_start_node(side_n_id_t i);
  void add_haplotype_stop_node(side_n_id_t i);

  void set_min_id(std::size_t min_id);
  void set_max_id(std::size_t max_id);

  // sort by in degree on the left side of the
  // void sort();

  std::map<std::size_t, component> count_components(const core::config& app_config);

  void set_vertex_eq_class(std::size_t v_idx, std::size_t eq_class);
  void set_edge_eq_class(std::size_t e_idx, std::size_t eq_class);

  // ----
  // misc
  // ----
  void dbg_print();
  void print_dot() const;


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

std::vector<VariationGraph> componetize(const VariationGraph &vg, const core::config& app_config);

struct component {
  id_t idx; // a unique identifier of the component

  std::set<id_t> starts; // tips at the start of a graph
  std::set<id_t> ends; // tips at the end of a graph

  // tips that are not haplotype start nodes. Will be deleted from the variation graph
  std::set<id_t> orphan_tips;

  std::set<id_t> non_tip_hap_starts; // haplotype start nodes that are not tips

  std::vector<id_t> vertices;

  VariationGraph vg;

  // ------------
  // Constructors
  // ------------

  component()
    : idx(0),
      starts(std::set<id_t>()), ends(std::set<id_t>()), orphan_tips(std::set<id_t>()), non_tip_hap_starts(std::set<id_t>()),
      vertices(std::vector<id_t>()) {}

  component(id_t idx)
    : idx(idx),
      starts(std::set<id_t>()), ends(std::set<id_t>()), orphan_tips(std::set<id_t>()), non_tip_hap_starts(std::set<id_t>()),
      vertices(std::vector<id_t>()) {}

  component(id_t idx,
            const std::set<id_t>& starts, const std::set<id_t>& ends, const std::set<id_t>& orphan_tips, const std::set<id_t>& non_tip_hap_starts,
            const std::vector<id_t>& vertices)
      : idx(idx), starts(starts), ends(ends), orphan_tips(orphan_tips), non_tip_hap_starts(non_tip_hap_starts),
        vertices(vertices) {}

  // -------
  // Getters
  // -------

  // friend declaration for operator<<
  friend std::ostream& operator<<(std::ostream& os, const component& comp);
};



}; // namespace bidirected
#endif
