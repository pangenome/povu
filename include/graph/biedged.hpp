#ifndef BIEDGED_HPP
#define BIEDGED_HPP
#include <cstddef>
#include <format>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <stack>


#include "../common/utils.hpp"
#include "../common/types.hpp"
#include "./bidirected.hpp"
#include "../common/types.hpp"
#include "./spanning_tree.hpp"


namespace biedged {
using namespace povu::graph_types;
namespace pc = povu::constants;
namespace pgt = povu::graph_types;
namespace pst = povu::spanning_tree;
namespace pbd = povu::bidirected;
namespace pu = povu::utils;
namespace pt = povu::types;
namespace pgt = povu::graph_types;


/**
 * (l , r) == (r, l)
 */
class Edge {
  std::size_t v1_idx;
  pgt::v_type_e v1_type;   // not necessary but should allow for easy methods
  std::size_t v2_idx;
  pgt::v_type_e v2_type;   // not necessary but should allow for easy methods
  pgt::color c;

  // TODO: remove?
  //std::string label; // or sequence only applicable to black edges
  std::size_t eq_class {pc::UNDEFINED_SIZE_T}; // equivalence class of the edge

public:
  // --------------
  // constructor(s)
  // --------------
  Edge();
  Edge(std::size_t v1_idx, VertexType v1_type, std::size_t v2_idx, VertexType v2_type, color c);


  // ---------
  // getter(s)
  // ---------
  std::size_t get_eq_class() const;
  color get_color() const;
  std::size_t get_v1_idx() const;
  std::size_t get_v2_idx() const;
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
  friend constexpr auto operator<=>(const Edge&, const Edge&) = default;
};

/**
   a vertex is
 */
class Vertex {
  // indexes to the edge vector in Graph
  std::size_t black_edge;
  std::set<std::size_t> grey_edges;

  pgt::v_type_e type; // is this a 5' or a 3' vertex

  /*
    the name of the vertex in the input GFA
    note that 2 BEC vertices have the same handle
   */
  //std::string handle; // rename to id and store a size_t
  std::size_t id_;

  // index of the vertex in the vertex vector
  std::size_t vertex_idx;

public:
  // --------------
  // constructor(s)
  // --------------
  Vertex();
  Vertex(std::size_t id, std::size_t vertex_idx, v_type_e vertex_type);


  // ---------
  // getter(s)
  // ---------
  bool is_reversed() const;
  std::size_t id() const;
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
   * * always adds dummy vertices
   * the dummy vertices are dummy start and stop nodes
   *
   * @g the bidirected variation graph
   * @add_dummy_vertices if false, do not dummy vertices to the graph
   */
  BVariationGraph(const pbd::VG& g);


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
  /**
    * @brief Add an edge to the graph
    * @param v1_idx index of the first vertex in the vertex vector
    * @param v1_type type of the first vertex
    * @param v2_idx index of the second vertex in the vertex vector
    * @param v2_type type of the second vertex
    * @param c color of the edge
    * @return the index of the edge in the edge vector
   */
  std::size_t add_edge(std::size_t v1_idx, v_type v1_type,
                       std::size_t v2_idx, v_type v2_type,
                       pgt::color_e c);

  void add_vertex(std::size_t v_id, VertexType vertex_type);
  void update_eq_classes(pst::Tree &tree);


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
  pst::Tree compute_spanning_tree() const;
};

typedef BVariationGraph BVG;

}; // namespace biedged
#endif
