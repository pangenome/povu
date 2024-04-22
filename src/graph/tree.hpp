#ifndef TREE_HPP
#define TREE_HPP

#include <set>
#include <string>
#include <vector>
#include <cstddef>
#include <map>
#include "../core/constants.hpp"
namespace tree {


/*
 * Vertex
 * ------
 *
 */
class Vertex {
  // TODO: make sure id is the same as the index in the di_graph or u_graph or have
  // a way of mapping to the index
  std::size_t id; // allow negative ids?
  std::size_t class_; // equivalence class

  std::size_t parent; // TODO: use underscore and make parent() and method
  std::set<std::size_t> children;

  // std::size_t depth; // depth of the vertex in the tree

  bool is_valid_; // is this a valid vertex in the tree
  bool is_dummy_node_; // is this a dummy node in the tree

  std::string meta;

public:
  // --------------
  // constructor(s)
  // --------------
  Vertex(); // a null vertex
  Vertex(std::size_t id); // root constructor sets the parent to a max value
  Vertex(std::size_t id, bool is_dummy);
  //Vertex(std::size_t id, std::size_t eq_class); // root constructor sets the parent to a max value
  Vertex(std::size_t id, std::size_t parent_id);

  // non-root vertex constructor
  Vertex(std::size_t id, std::size_t parent_id, std::size_t eq_class);

  Vertex(std::size_t id, std::size_t parent_id, std::size_t eq_class, bool is_dummy);

  // ---------
  // getter(s)
  // ---------
  bool is_valid() const;
  bool is_dummy() const;
  std::set<std::size_t> const& get_children() const;
  std::size_t get_id() const;
  std::size_t get_class() const;
  std::size_t get_parent() const;
  std::string get_meta() const;


 // ---------
  // setter(s)
  // ---------
  void set_parent(std::size_t child_id);
  void add_child(std::size_t child_id);
  void remove_child(std::size_t child_id);
  void set_class(std::size_t class_);
  void set_meta(std::string&& meta);
};

/*
 * Tree
 * ------
 *
 * a simple vector backed tree implementation
 */
class Tree {
  std::vector<Vertex> vertices;
  //static const std::size_t root_idx_{}; // root idx can change
  std::size_t root_idx_;
public:
  // --------------
  // constructor(s)
  // --------------
  // construct a null tree
  // without any vertices
  Tree();
  // construct a tree with n null vertices but a valid root
  Tree(std::size_t n, bool artificial_root=false);


  // ---------
  // getter(s)
  // ---------
  Vertex const& get_root() const;

  // given the id of the vertex
  // returns the set of children of the vertex
  std::set<std::size_t> const& get_children(std::size_t id) const;

  std::size_t get_class(std::size_t id) const;
  std::size_t get_meta(std::size_t id) const;
  Vertex const& get_vertex(std::size_t id) const;
  Vertex& get_vertex_mut(std::size_t id);
  std::size_t get_parent(std::size_t id) const;

  // the number of vertices in the tree
  std::size_t size() const;

  bool empty() const;


  // ---------
  // setter(s)
  // ---------

  // add a vertex to the tree
  // given the id of the vertex and the id of the parent
  // returns false if the vertex is already in the tree or some error occurs
  // returns true if the vertex is added to the tree
  // TODO: deprecated
  bool add_vertex(std::size_t parent_id, std::size_t id);
  bool add_vertex(std::size_t parent_id, std::size_t id, std::size_t eq_class);
  bool add_vertex(std::size_t parent_id, std::size_t id, std::size_t eq_class, std::string& meta);
  bool add_vertex(std::size_t parent_id, std::size_t id, std::size_t eq_class, std::string& meta, bool is_dummy);

  // add a root to the tree
  //
  // bool add_root(std::size_t id, std::size_t eq_class, std::string& meta, bool is_dummy);

  bool remove_vertex(std::size_t id);

  bool set_root(std::size_t id);

  // -----------------
  // display method(s)
  // -----------------

  // dot format output of the tree
  void print_dot(bool with_classes=false);
};

struct node {
    std::set<std::size_t> children;
    std::size_t parent {core::constants::UNDEFINED_SIZE_T};
};


/**
 * Map backed tree
 * ---------------
 *
 */
template <typename KeyType, typename ValueType>
using MapTree = std::map<KeyType, ValueType>;

} // namespace tree

#endif
