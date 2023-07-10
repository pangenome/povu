#ifndef VST_HPP
#define VST_HPP

#include "../graph/spanning_tree.hpp"
#include <cstddef>

namespace vst {

void cycle_equiv(spanning_tree::Tree &t);

class Vertex {
  std::size_t idx;
  std::size_t cycle_equiv_class_;

  std::size_t parent;
  std::set<std::size_t> children;

public:

  // constructor(s)
  // --------------
  Vertex(std::size_t idx); // root
  
  // root for use by PST
  Vertex(std::size_t idx, std::size_t cycle_equiv_class_); 
  Vertex(std::size_t idx, std::size_t cycle_equiv_class_, std::size_t parent);

  // getters
  // -------
  // given the idx of the child in the tree vector
  std::set<std::size_t> const &get_children() const;

  // setters
  // -------
  // given the idx of the child in the tree vector
  void add_child(std::size_t child_idx);
};

class VST {
  std::vector<Vertex> t;
  static const std::size_t root_idx_{};

public:
  // constructor(s)
  // --------------
  VST(spanning_tree::Tree &t);

  // getters
  // -------
  std::size_t size() const;
  std::size_t root_idx() const;
  std::set<std::size_t> const& get_children(std::size_t v) const;

  // IO
  // --
  void print_dot();
};
} // namespace vst


namespace tree {

// Vertex
// ======
  
class Vertex {
  std::size_t id;
  std::size_t parent;
  std::set<std::size_t> children;

  bool is_valid_; // is this a valid vertex in the tree

public:
  // constructor(s)
  // ----------------
  Vertex(); // a null vertex
  Vertex(std::size_t id); // root constructor sets the parent to a max value
  Vertex(std::size_t id, std::size_t parent_id); // non-root vertex constructor

    
  // getters
  // -------
  bool is_valid() const;
  std::set<std::size_t> const& get_children() const;

  // setters
  // -------
  void add_child(std::size_t child_id);
};

// Tree
// ====
  
// a simple vector backed tree implementation 
class Tree {
  std::vector<Vertex> vertices;
  static const std::size_t root_idx_{}; // root is always at index 0

public:
  // constructor(s)
  // --------------
  Tree(); // construct a null tree
  Tree(std::size_t n); // construct a tree with n null vertices but a valid root

  // setters
  // -------

  // add a vertex to the tree
  // given the id of the vertex and the id of the parent
  // returns false if the vertex is already in the tree or some error occurs
  // returns true if the vertex is added to the tree
  bool add_vertex(std::size_t parent_id, std::size_t id);

  // getters
  // -------

  // given the id of the vertex
  // returns the set of children of the vertex
  std::set<std::size_t> const& get_children(std::size_t id) const;

  // the number of vertices in the tree
  std::size_t size() const;

  // dot format output of the tree
  void print_dot();
};

} // namespace tree


namespace pst {
// a PST as defined in the paper
// -----------------------------

  

class PST {
  std::vector<vst::Vertex> t;
  static const std::size_t root_idx_{};

public:
  // constructor(s)
  // --------------
  // given a spanning tree with the equivalence classes generate the PST
  PST(spanning_tree::Tree &t);

  // getters
  // -------
  std::size_t size() const;
  std::size_t root_idx() const;
  std::set<std::size_t> const& get_children(std::size_t v) const;

  // IO
  // --
  void print_dot();
};

  tree::Tree compute_pst(spanning_tree::Tree &st);
  
} // namespace pst


#endif
