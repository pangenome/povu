#include <iostream>
#include <format>
#include <limits>
#include "./tree.hpp"

namespace tree {
  // TODO: declare all of these in a header file
  const std::size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();

  // Vertex
  // ======

  // constructor(s)
  // --------------
  Vertex::Vertex() :
    id(SIZE_T_MAX), parent(SIZE_T_MAX), children(std::set<std::size_t>{}), is_valid_(false) {}
  Vertex::Vertex(std::size_t id) :
    id(id), parent(SIZE_T_MAX), children(std::set<std::size_t>{}), is_valid_(true) {}
  Vertex::Vertex(std::size_t id, std::size_t parent_id) :
    id(id), parent(parent_id), children(std::set<std::size_t>{}), is_valid_(true) {}

  // member function(s)
  // ------------------

  // getters
  // -------
  bool Vertex::is_valid() const { return this->is_valid_; }
  std::set<std::size_t> const& Vertex::get_children() const { return this->children; }

  // setters
  // -------

  void Vertex::add_child(std::size_t child_id) {
    this->children.insert(child_id);
  }

  // Tree
  // ====

  // constructor(s)

  Tree::Tree() : vertices(std::vector<tree::Vertex>{}) {}

  Tree::Tree(std::size_t n) : vertices(std::vector<Vertex>(n, Vertex())) {}

  // member function(s)
  // ------------------

  // setters
  bool Tree::add_vertex(std::size_t parent_id, std::size_t id) {
    // TODO: there's a logical error in the caller if the vertex is already in the tree
    //       should we throw an exception here?
    if (this->vertices[id].is_valid()) { return false;  }
    this->vertices[id] = Vertex(id, parent_id);
    this->vertices[parent_id].add_child(id);
    return true;
  }

  std::set<std::size_t> const& Tree::get_children(std::size_t v) const {
    return this->vertices.at(v).get_children();
  }

  void Tree::print_dot() {
  std::cout << std::format(
    "graph G {{\n"
    "\trankdir = TB;\n"
    "\tnode[shape = circle];\n"
    "\tedge [arrowhead=vee];\n"
  );

  for (std::size_t i{}; i < this->size(); i++) {
    for (auto c : this->get_children(i)) {
      std::cout << std::format("\t{} -- {};\n", i, c);
    }
  }

  std::cout << "}" << std::endl;
}

  // getters

  // TODO: should this be the number of valid vertices?
  std::size_t Tree::size() const { return this->vertices.size(); }

} // namespace tree
