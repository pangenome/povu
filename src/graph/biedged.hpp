#ifndef BIEDGED_HPP
#define BIEDGED_HPP

#include "../core/core.hpp"
#include "./bidirected.hpp"

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
	r
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

  std::string label; // or sequence only applicable to black edges

public:
  
  // ------------
  // constructors
  // ------------
  
  Edge();
  Edge(std::size_t v1, VertexType v1_type, std::size_t v2, VertexType v2_type, core::color c);
  Edge(std::size_t v1, VertexType v1_type, std::size_t v2, VertexType v2_type, core::color c, std::string label);


  
// spaceship operator
friend constexpr auto operator<=>(Edge, Edge) = default;
  
  
  // -------
  // getters
  // -------
  const std::string& get_label() const;
  std::size_t get_v1_idx() const;
  std::size_t get_v2_idx() const;
  core::color get_color() const;

  // implement << operator
  friend std::ostream& operator<<(std::ostream& os, const Edge& e);
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

  // index of the vertex in the vertex vector
  std::size_t vertex_idx;
  
public:

  // ------------
  // constructors
  // ------------
  
  Vertex();
  Vertex(const std::string& id, std::size_t vertex_idx, VertexType vertex_type);

  /*
	setters and getters
	-------------------
  */

  // -------
  // getters
  // -------
  
  bool is_reversed() const;
  const std::string& get_handle() const;
  VertexType get_type() const;

  // -------
  // setters
  // -------
  
  // returns the new value
  bool toggle_reversed();
  void add_edge(std::size_t edge_index, core::color c);

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
  BVariationGraph(const bidirected::VariationGraph& g);

  // -------
  // getters
  // -------
  Vertex& get_vertex_mut(std::size_t vertex_idx);
  const Vertex& get_vertex(std::size_t vertex_idx) const;

  // -------
  // setters
  // -------
    
  void add_edge(std::size_t v1,  VertexType v1_type,
				std::size_t v2, VertexType v2_type,
				core::color c);
  
  void add_edge(std::size_t v1, VertexType v1_type,
				std::size_t v2, VertexType v2_type,
				core::color c, std::string label);

  void add_vertex(std::string handle_str, VertexType vertex_type);
  bool replace_vertex(std::size_t vertex_idx, std::string handle_str, VertexType vertex_type);

  // ---------------
  // display methods
  // ---------------
  
  void print_dot() const;
};
  
}; // namespace biedged
#endif
