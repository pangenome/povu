#ifndef TREE_HPP
#define TREE_HPP

#include <cstddef>
#include <format>
#include <iostream>
#include <string>
#include <sys/types.h>
#include <tuple>
#include <unistd.h>
#include <utility>
#include <vector>
#include <variant>
#include <optional>

#include "../common/types.hpp"

// generic tree implementation
namespace povu::tree {
using namespace povu::constants;


template <typename T> class Vertex {
  std::size_t id;
  std::optional<T> data_;

public:
  // --------------
  // constructor(s)
  // --------------
  Vertex(std::size_t id) : id(id), data_(std::nullopt) {}
  Vertex(std::size_t id, T d) : id(id), data_(d) {}
  //Vertex(std::size_t id, id_n_cls r) : id(id), data_({r}) {}

  // ---------
  // getter(s)
  // ---------
  std::size_t get_id() const { return id; }
  //  VertexType get_type() const { return type_; }
  std::optional<T> get_data() const { return data_; }

  // TODO: remove this?
  std::string data_as_str() const {
    if (data_.has_value()) {
      return "..";
    }
    else {
      return ".";
    }
  }
  std::string as_str() const {
    return std::format("v{}:{}", this->get_id(), this->data_as_str());
  }
  // ---------
  // setter(s)
  // ---------
  void update_data(T new_data) {
    this->data_ = new_data;
  }
};


// always has a dummy root vertex
template <typename T>  class Tree {
  std::vector<Vertex<T>> vertices;
  std::vector<std::size_t> parent_v; // parent of each vertex
  std::vector<std::vector<std::size_t>> children_v; // children of each vertex
  pt::idx_t root_idx_; // index of the root vertex in the vertices vector

public:
  // --------------
  // constructor(s)
  // --------------

  Tree() {
    vertices.push_back(Vertex<T>(INVALID_ID));
    root_idx_ = vertices.size() - 1;
  }

  Tree(std::size_t expected_size) : Tree() {
    vertices.reserve(expected_size);
    this->parent_v = std::vector<std::size_t>(expected_size+1, INVALID_ID);
    this->children_v = std::vector<std::vector<std::size_t>>(1+expected_size);
  }

  // ---------
  // getter(s)
  // ---------
  // deprecated
  [[deprecated("use vtx_count() instead")]]
  std::size_t size() const {
    return this->vertices.size();
  }

  std::size_t vtx_count() const {
    return this->vertices.size();
  }

  pt::idx_t root_idx() const {
    return this->root_idx_;
  }

  const Vertex<T>& get_root() const {
    return this->vertices[this->root_idx()];
  }

  const Vertex<T>& get_vertex(std::size_t v_idx) const {
    return this->vertices[v_idx];
  }

  Vertex<T>& get_vertex_mut(std::size_t v_idx) {
    return this->vertices[v_idx];
  }

  const Vertex<T>& get_parent(std::size_t v_idx) const {
    return this->vertices[parent_v[v_idx]];
  }

  std::size_t get_parent_idx(std::size_t v_idx) const {
    return this->parent_v[v_idx];
  }

  bool is_leaf(std::size_t v_idx) const {
    return  v_idx >= this->children_v.size() || this->children_v[v_idx].empty();
  }

  const std::vector<std::size_t>& get_children(std::size_t v_idx) const {
    return this->children_v[v_idx];
  }

  // ---------
  // setter(s)
  // ---------

  /**
    * @brief Add a vertex to the tree
    * @param v Vertex to be added
    * @return Index of the added vertex
   */
  std::size_t add_vertex(Vertex<T> v) {
    this->vertices.push_back(v);
    return this->vertices.size() - 1;
  }

  void add_edge(std::size_t parent, std::size_t child) {

    while (child >= this->parent_v.size()) { this->parent_v.push_back(INVALID_ID); }
    this->parent_v[child] = parent;
    while (parent >= this->children_v.size()) { this->children_v.push_back(std::vector<std::size_t>()); }
    this->children_v[parent].push_back(child);
  }

  // ----
  // misc
  // ----
  void print_dot() {
  std::cout << std::format(
    "graph G {{\n"
    "\trankdir = TD;\n"
    "\tnode[shape = circle];\n"
    "\tedge [arrowhead=vee];\n"
  );

  // print vertices
  for (std::size_t i{}; i < this->size(); i++) {
  }

  // print edges
  for (std::size_t i{}; i < this->size(); i++) {
    for (std::size_t c : this->get_children(i)) {
      std::cout << std::format("\t{} -- {};\n", i, c);
    }
  }

  std::cout << "}" << std::endl;
}

};


} // namespace tree




#endif
