#ifndef BIEDGED_HPP
#define BIEDGED_HPP

//#include "../common/utils.hpp"
#include "./bidirected.hpp"
#include "./spanning_tree.hpp"
#include "../common/common.hpp"
#include <cstddef>


namespace biedged {
using namespace graph_types;
namespace constants = common::constants;

/**
   (l , r) == (r, l)
 */
class Edge {
  std::size_t v1_idx;
  VertexType v1_type;   // not necessary but should allow for easy methods
  std::size_t v2_idx;
  VertexType v2_type;   // not necessary but should allow for easy methods
  color c;

  // TODO: remove?
  std::string label; // or sequence only applicable to black edges
  std::size_t eq_class {constants::UNDEFINED_SIZE_T}; // equivalence class of the edge

public:
  // --------------
  // constructor(s)
  // --------------
  Edge();
  Edge(std::size_t v1, VertexType v1_type, std::size_t v2, VertexType v2_type, color c);
  Edge(std::size_t v1, VertexType v1_type, std::size_t v2, VertexType v2_type, color c, std::string label);


  // ---------
  // getter(s)
  // ---------
  std::size_t get_eq_class() const;
  color get_color() const;
  std::size_t get_v1_idx() const;
  std::size_t get_v2_idx() const;
  const std::string& get_label() const;
  std::size_t get_other_vertex(std::size_t vertex_index) const;


  // ---------
  // setter(s)
  // ---------
  void set_eq_class(std::size_t e);
  void set_v1_idx(std::size_t i);
  void set_v2_idx(std::size_t i);

  // -----------
  // operator(s)
  // -----------
  // implement << operator
  friend std::ostream& operator<<(std::ostream& os, const Edge& e);
  // spaceship operator
  friend constexpr auto operator<=>(Edge, Edge) = default;
};

/**
   a vertex is
 */
class Vertex {
  // indexes to the edge vector in Graph
  std::size_t black_edge;
  std::set<std::size_t> grey_edges;

  VertexType type; // is this a 5' or a 3' vertex

  /*
    the name of the vertex in the input GFA
    note that 2 BEC vertices have the same handle
   */
  std::string handle;

  // index of the vertex in the vertex vector
  std::size_t vertex_idx;

public:
  // --------------
  // constructor(s)
  // --------------
  Vertex();
  Vertex(const std::string& id, std::size_t vertex_idx, VertexType vertex_type);


  // ---------
  // getter(s)
  // ---------
  bool is_reversed() const;
  const std::string& get_handle() const;
  VertexType get_type() const;
  std::set<std::size_t> const& get_grey_edges() const;
  std::size_t get_black_edge() const;
  std::size_t get_vertex_idx() const;
  std::vector<std::size_t> get_neighbours() const;


  // ---------
  // setter(s)
  // ---------
  // returns the new value
  bool toggle_reversed();
  void add_edge(std::size_t edge_index, color c);
  void set_vertex_idx(std::size_t i);

  // It is up to the user to make sure that the path_id is not already in the "set"
  //void add_path(std::size_t path_id, std::size_t step_index);
};


class BVariationGraph {
  std::vector<Vertex> vertices;
  std::vector<Edge> edges;


  // indexes of the dummy vertices in the vertex vector
  std::vector<std::size_t> dummy_vertices_;


public:
  // ------------
  // constructors
  // ------------
  /**
    * @brief Construct a bi-edged graph from a bidirected variation graph
    *
    * the dummy vertices are dummy start and stop nodes
    *
    *
    *
    *
    * @g the bidirected variation graph
    * @add_dummy_vertices if false, do not dummy vertices to the graph
    */
  BVariationGraph(const bidirected::VariationGraph& g, bool add_dummy_vertices=true);


  // -------
  // getters
  // -------
  std::size_t size() const;
  std::size_t num_edges() const;

  Edge &get_edge_mut(std::size_t edge_idx);
  const Edge& get_edge(std::size_t edge_idx) const;
  const std::vector<Edge>& get_all_edges() const;
  Vertex& get_vertex_mut(std::size_t vertex_idx);
  const Vertex& get_vertex(std::size_t vertex_idx) const;
  const std::vector<size_t> &get_dummy_vertices() const;
  std::vector<std::pair<color, std::size_t>> get_neighbours(std::size_t vertex_idx) const;


  // -------
  // setters
  // -------
  std::size_t add_edge(std::size_t v1, VertexType v1_type,
                       std::size_t v2, VertexType v2_type,
                       color c);

  std::size_t add_edge(std::size_t v1, VertexType v1_type,
                       std::size_t v2, VertexType v2_type,
                       color c, std::string label);

  void add_vertex(std::string handle_str, VertexType vertex_type);
  bool replace_vertex(std::size_t vertex_idx, std::string handle_str, VertexType vertex_type);
  void update_eq_classes(spanning_tree::Tree &tree);


  // ---------------
  // display methods
  // ---------------
  void print_dot() const;


  // ----
  // misc
  // ----

  /**
   * @brief Compute the DFS spanning tree of the graph
   */
  spanning_tree::Tree compute_spanning_tree() const;
};

}; // namespace biedged
#endif
