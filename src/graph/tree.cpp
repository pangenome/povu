#include <iostream>
#include <format>
#include <limits>
#include <string>

#include "./tree.hpp"
#include "../core/core.hpp"

namespace tree {
  // TODO: use the one in constants.hpp
  const std::size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();

  // Vertex
  // ======

  // constructor(s)
  // --------------
  Vertex::Vertex() :
    id(core::constants::UNDEFINED_SIZE_T),
    class_(core::constants::UNDEFINED_SIZE_T),
    parent(core::constants::UNDEFINED_SIZE_T),
    children(std::set<std::size_t>{}),
    is_valid_(false) {}

  Vertex::Vertex(std::size_t id) :
    id(id),
    class_(core::constants::UNDEFINED_SIZE_T),
    parent(core::constants::UNDEFINED_SIZE_T),
    children(std::set<std::size_t>{}),
    is_valid_(true) {}

  Vertex::Vertex(std::size_t id, std::size_t parent_id) :
    id(id),
    class_(core::constants::UNDEFINED_SIZE_T),
    parent(parent_id),
    children(std::set<std::size_t>{}),
    is_valid_(true) {}

  Vertex::Vertex(std::size_t id, std::size_t parent_id, std::size_t eq_class) :
    id(id),
    class_(eq_class),
    parent(parent_id),
    children(std::set<std::size_t>{}),
    is_valid_(true) {}

  // member function(s)
  // ------------------

  // getters
  // -------
  bool Vertex::is_valid() const { return this->is_valid_; }
  std::set<std::size_t> const& Vertex::get_children() const { return this->children; }
std::size_t Vertex::get_class() const { return this->class_; }
std::size_t Vertex::get_parent() const { return this->parent; }
  // setters
  // -------

  void Vertex::add_child(std::size_t child_id) {
    this->children.insert(child_id);
  }


  void Vertex::remove_child(std::size_t child_id) {
    this->children.erase(child_id);
  }
  
  // Tree
  // ====

  // constructor(s)

  Tree::Tree() : vertices(std::vector<tree::Vertex>{}) {}

  Tree::Tree(std::size_t n, bool artificial_root) : vertices(std::vector<Vertex>(n, Vertex())) {
    if (artificial_root) {
      this->vertices[0] = Vertex(0);
    }
  }

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

bool Tree::add_vertex(std::size_t parent_id, std::size_t id, std::size_t eq_class) {
    // TODO: there's a logical error in the caller if the vertex is already in the tree
    //       should we throw an exception here?
    if (this->vertices[id].is_valid()) { return false;  }
    this->vertices[id] = Vertex(id, parent_id, eq_class);
    this->vertices[parent_id].add_child(id);
    return true;
}

bool Tree::remove_vertex(std::size_t id) {
    if (!this->vertices[id].is_valid()) { return false; }
    std::size_t parent_id = this->vertices[id].get_parent();
    this->vertices[parent_id].remove_child(id);
    this->vertices[id] = Vertex();
    return true;
}
  
std::set<std::size_t> const& Tree::get_children(std::size_t v) const {
  return this->vertices.at(v).get_children();
}

std::size_t Tree::get_parent(std::size_t id) const {
  std::size_t p= this->vertices.at(id).get_parent();
  return p == core::constants::UNDEFINED_SIZE_T ? 0 : p; 
  //return this->vertices.at(id).get_parent();
}

Vertex const& Tree::get_vertex(std::size_t id) const {
  
  return this->vertices.at(id);
}
  
void Tree::print_dot(bool with_classes) {
  std::cout << std::format(
    "graph G {{\n"
    "\trankdir = TB;\n"
    "\tnode[shape = circle];\n"
    "\tedge [arrowhead=vee];\n"
  );


  if (with_classes) {
    for (std::size_t i{}; i < this->size(); i++) {

      // TODO: remove
      if (this->vertices[i].get_class() == core::constants::UNDEFINED_SIZE_T) { continue; }
      
      std::string class_label =
        this->vertices[i].get_class() == core::constants::UNDEFINED_SIZE_T ?
        "UNDEFINED" : std::to_string(this->vertices[i].get_class());

      
      
      std::cout << std::format("\t{} [label=\"v: {}\\ncl: {}\"];\n", i, i, class_label);
    }
  }

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
