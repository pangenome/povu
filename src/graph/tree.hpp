#ifndef TREE_HPP
#define TREE_HPP


#include <cstddef>
#include <set>
#include <vector>


namespace tree {

// Vertex
// ======

class Vertex {
  std::size_t id; // allow negative ids?
  std::size_t class_; // equivalence class

  std::size_t parent;
  std::set<std::size_t> children;

  // std::size_t depth; // depth of the vertex in the tree

  bool is_valid_; // is this a valid vertex in the tree

public:
  // constructor(s)
  // --------------
  Vertex(); // a null vertex
  Vertex(std::size_t id); // root constructor sets the parent to a max value
  Vertex(std::size_t id, std::size_t parent_id);

  // non-root vertex constructor
  Vertex(std::size_t id, std::size_t parent_id, std::size_t eq_class); // non-root vertex constructor

  // getters
  // -------
  bool is_valid() const;
  std::set<std::size_t> const& get_children() const;
  std::size_t get_class() const;

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
  // construct a null tree
  Tree();
  // construct a tree with n null vertices but a valid root
  Tree(std::size_t n, bool artificial_root=false);

  // setters
  // -------

  // add a vertex to the tree
  // given the id of the vertex and the id of the parent
  // returns false if the vertex is already in the tree or some error occurs
  // returns true if the vertex is added to the tree
  // TODO: deprecated
  bool add_vertex(std::size_t parent_id, std::size_t id);

  bool add_vertex(std::size_t parent_id, std::size_t id, std::size_t eq_class);

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

#endif
