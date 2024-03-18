#ifndef BIEDGED_HPP
#define BIEDGED_HPP

#include "../core/core.hpp"
#include "./bidirected.hpp"
#include "./spanning_tree.hpp"
#include <cstddef>

namespace biedged {

// Merge path_t and PathInfo into one namespace
/**
 * Path (or color)
 * ---------------
 */
struct path_t {
  std::string name; // name as pertains the GFA file
  std::size_t id; // numerical id to be associated with handle
  bool is_circular; // is the path circular?
};


struct unordered_pair{
  std::size_t l;
  std::size_t r;

  unordered_pair(std::size_t l,std::size_t r):
    l(std::min(l,r)), r(std::max(l,r)) {}

  // spaceship operator
  friend constexpr auto operator<=>(unordered_pair, unordered_pair) = default;
};

struct PathInfo {
  std::size_t path_id;
  std::size_t step_index;

  PathInfo(): path_id(0), step_index(0) {}
  PathInfo(std::size_t path_id, std::size_t step_index): path_id(path_id), step_index(step_index) {}
};


/**
   l (left) 5' or +
   r (left) 3' or -
 */
enum class VertexType {
    l,
    r,
    dummy
};

// implement << operator for VertexType
std::ostream& operator<<(std::ostream& os, const VertexType& vt);

/**
   (l , r) == (r, l)
 */
class Edge {
  std::size_t v1_idx;
  VertexType v1_type;   // not necessary but should allow for easy methods
  std::size_t v2_idx;
  VertexType v2_type;   // not necessary but should allow for easy methods
  core::color c;

  // TODO: remove?
  std::string label; // or sequence only applicable to black edges

public:
  // ------------
  // constructors
  // ------------
  Edge();
  Edge(std::size_t v1, VertexType v1_type, std::size_t v2, VertexType v2_type, core::color c);
  Edge(std::size_t v1, VertexType v1_type, std::size_t v2, VertexType v2_type, core::color c, std::string label);


  // -------
  // getters
  // -------
  const std::string& get_label() const;
  std::size_t get_v1_idx() const;
  std::size_t get_v2_idx() const;
  core::color get_color() const;


  // -------
  // setters
  // -------
  void set_v1_idx(std::size_t i);
  void set_v2_idx(std::size_t i);

  // ---------
  // operators
  // ---------
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

  // paths (also colors)
  // TODO: change to unordered_set
  std::vector<PathInfo> paths;

  VertexType type; // is this a 5' or a 3' vertex

  // from libHandleGraph
  // 2 BEC vertices have the same handle
  std::string handle;
  bool is_reversed_;
  // std::size_t bidirected_idx // index of the vertex in the bidirected graph

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
  std::set<std::size_t> get_grey_edges() const;
  std::size_t get_black_edge() const;
  std::size_t get_vertex_idx() const;


  // ---------
  // setter(s)
  // ---------
  // returns the new value
  bool toggle_reversed();
  void add_edge(std::size_t edge_index, core::color c);
  void set_vertex_idx(std::size_t i);

  // It is up to the user to make sure that the path_id is not already in the "set"
  void add_path(std::size_t path_id, std::size_t step_index);
};


class BVariationGraph {
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
  Vertex& get_vertex_mut(std::size_t vertex_idx);
  const Vertex& get_vertex(std::size_t vertex_idx) const;
  std::size_t size() const;


  // -------
  // setters
  // -------
  std::size_t add_edge(std::size_t v1,  VertexType v1_type,
                       std::size_t v2, VertexType v2_type,
                       core::color c);

  std::size_t add_edge(std::size_t v1, VertexType v1_type,
                std::size_t v2, VertexType v2_type,
                core::color c, std::string label);

  void add_vertex(std::string handle_str, VertexType vertex_type);
  bool replace_vertex(std::size_t vertex_idx, std::string handle_str, VertexType vertex_type);


  // ---------------
  // display methods
  // ---------------
  void print_dot() const;


  // ----
  // misc
  // ----

  // If the graph is not a single SCC, then make it one.
  void componetize(); // TODO: remove
  spanning_tree::Tree compute_spanning_tree() const;
};

}; // namespace biedged
#endif
