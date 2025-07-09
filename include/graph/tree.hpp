#ifndef TREE_HPP
#define TREE_HPP

#include <cstddef>
#include <format>
#include <iostream>
#include <optional>
#include <string>
#include <sys/types.h>
#include <tuple>
#include <type_traits>
#include <unistd.h>
#include <utility>
#include <variant>
#include <vector>

#include "../common/types/types.hpp"

// generic tree implementation
namespace povu::tree {
using namespace povu::constants;
namespace pvst = povu::types::pvst;

// always has a dummy root vertex
class Tree {
  std::vector<std::unique_ptr<pvst::VertexBase>> vertices;
  std::vector<pt::idx_t> parent_v; // parent of each vertex
  std::vector<std::vector<pt::idx_t>> children_v; // children of each vertex
  pt::idx_t root_idx_; // index of the root vertex in the vertices vector

public:
  // --------------
  // constructor(s)
  // --------------

  Tree(): root_idx_(INVALID_IDX) {}

  Tree(pt::idx_t expected_size) : Tree() {
    vertices.reserve(expected_size);
    this->parent_v = std::vector<pt::idx_t>(expected_size+1, INVALID_ID);
    this->children_v = std::vector<std::vector<pt::idx_t>>(1+expected_size);
  }

  // ---------
  // getter(s)
  // ---------

  // [[deprecated("use vtx_count() instead")]]
  // pt::idx_t size() const {
  //   return this->vertices.size();
  // }

  pt::idx_t vtx_count() const {
    return this->vertices.size();
  }

  pt::idx_t root_idx() const {
    return this->root_idx_;
  }

  const pvst::VertexBase& get_root() const {
    return this->get_vertex(this->root_idx());
  }

  const pvst::VertexBase &get_vertex(pt::idx_t v_idx) const {
    return *this->vertices[v_idx];
  }

  pvst::VertexBase& get_vertex_mut(pt::idx_t v_idx) {
    return *this->vertices[v_idx];
  }

  const pvst::VertexBase& get_parent(pt::idx_t v_idx) const {
    return *this->vertices[parent_v[v_idx]];
  }

  pt::idx_t get_parent_idx(pt::idx_t v_idx) const {
    return this->parent_v[v_idx];
  }

  bool is_leaf(pt::idx_t v_idx) const {
    return  v_idx >= this->children_v.size() || this->children_v[v_idx].empty();
  }

  const std::vector<pt::idx_t>& get_children(pt::idx_t v_idx) const {
    return this->children_v[v_idx];
  }



  // ---------
  // setter(s)
  // ---------

  void set_root_idx(pt::idx_t v_idx) {
    if (this->root_idx_ != INVALID_IDX){
      throw std::logic_error("Root index is already set");
    }

    if (v_idx >= this->vertices.size()) {
      throw std::out_of_range("Vertex index out of range");
    }

    this->root_idx_ = v_idx;
  }

  /**
    * @brief Add a vertex to the tree
    * @param v Vertex to be added
    * @return Index of the added vertex
   */
  template <typename T> pt::idx_t add_vertex(T v) {
    pt::idx_t v_idx = this->vertices.size();
    v.set_idx(v_idx); // set the index of the vertex

    while (v_idx >= this->parent_v.size()) {
      this->parent_v.push_back(INVALID_ID); // ensure parent_v has enough space
    }

    while (v_idx >= this->children_v.size()) {
      this->children_v.push_back(std::vector<pt::idx_t>()); // ensure children_v has enough space
    }

    // create a unique pointer to the vertex and add it to the vertices vector
    auto ptr = std::make_unique<T>(v);
    this->vertices.push_back(std::move(ptr));
    
    return v_idx;
  }

  void add_edge(pt::idx_t parent, pt::idx_t child) {
    while (child >= this->parent_v.size()) { this->parent_v.push_back(INVALID_ID); }
    this->parent_v[child] = parent;
    while (parent >= this->children_v.size()) { this->children_v.push_back(std::vector<pt::idx_t>()); }
    this->children_v[parent].push_back(child);
  }

  void del_edge(pt::idx_t parent, pt::idx_t child) {
    this->parent_v[child] = INVALID_IDX;

    std::vector<pt::idx_t>&children = this->children_v[parent];
    auto it = std::find(children.begin(), children.end(), child);
    if (it != children.end()) {
      children.erase(it);
    }
  }

  // ----
  // misc
  // ----
  void print_dot() const {
  std::cout << std::format(
    "graph G {{\n"
    "\trankdir = TD;\n"
    "\tnode[shape = circle];\n"
    "\tedge [arrowhead=vee];\n"
  );

  // print vertices
  for (pt::idx_t i{}; i < this->vtx_count(); i++) {
    std::cout << std::format("\t{} [label=\"{}\"];\n", i, this->get_vertex(i).as_str());
  }

  // print edges
  for (pt::idx_t i{}; i < this->vtx_count(); i++) {
    for (pt::idx_t c : this->get_children(i)) {
      std::cout << std::format("\t{} -- {};\n", i, c);
    }
  }

  std::cout << "}" << std::endl;
}

};


} // namespace povu::tree

namespace povu::tree_old {
using namespace povu::constants;


template <typename T> class Vertex {
  pt::idx_t id;
  std::optional<T> data_;

public:
  // --------------
  // constructor(s)
  // --------------
  Vertex(pt::idx_t id) : id(id), data_(std::nullopt) {}
  Vertex(pt::idx_t id, T d) : id(id), data_(d) {}
  //Vertex(std::size_t id, id_n_cls r) : id(id), data_({r}) {}

  // ---------
  // getter(s)
  // ---------
  pt::idx_t get_id() const { return id; }
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
  std::vector<T> vertices;
  std::vector<pt::idx_t> parent_v; // parent of each vertex
  std::vector<std::vector<pt::idx_t>> children_v; // children of each vertex
  pt::idx_t root_idx_; // index of the root vertex in the vertices vector

public:
  // --------------
  // constructor(s)
  // --------------

  Tree(): root_idx_(INVALID_IDX) {}

  Tree(pt::idx_t expected_size) : Tree() {
    vertices.reserve(expected_size);
    this->parent_v = std::vector<pt::idx_t>(expected_size+1, INVALID_ID);
    this->children_v = std::vector<std::vector<pt::idx_t>>(1+expected_size);
  }

  // ---------
  // getter(s)
  // ---------

  [[deprecated("use vtx_count() instead")]]
  pt::idx_t size() const {
    return this->vertices.size();
  }

  pt::idx_t vtx_count() const {
    return this->vertices.size();
  }

  pt::idx_t root_idx() const {
    return this->root_idx_;
  }

  const T& get_root() const {
    return this->get_vertex(this->root_idx());
  }

  const T& get_vertex(pt::idx_t v_idx) const {
    return this->vertices[v_idx];
  }

  T& get_vertex_mut(pt::idx_t v_idx) {
    return this->vertices[v_idx];
  }

  const T& get_parent(pt::idx_t v_idx) const {
    return this->vertices[parent_v[v_idx]];
  }

  pt::idx_t get_parent_idx(pt::idx_t v_idx) const {
    return this->parent_v[v_idx];
  }

  bool is_leaf(pt::idx_t v_idx) const {
    return  v_idx >= this->children_v.size() || this->children_v[v_idx].empty();
  }

  const std::vector<pt::idx_t>& get_children(pt::idx_t v_idx) const {
    return this->children_v[v_idx];
  }



  // ---------
  // setter(s)
  // ---------

  void set_root_idx(pt::idx_t v_idx) {
    if (this->root_idx_ != INVALID_IDX){
      throw std::logic_error("Root index is already set");
    }

    if (v_idx >= this->vertices.size()) {
      throw std::out_of_range("Vertex index out of range");
    }

    this->root_idx_ = v_idx;
  }

  /**
    * @brief Add a vertex to the tree
    * @param v Vertex to be added
    * @return Index of the added vertex
   */
  pt::idx_t add_vertex(T v) {
    pt::idx_t v_idx = this->vertices.size();

    while (v_idx >= this->parent_v.size()) {
      this->parent_v.push_back(INVALID_ID); // ensure parent_v has enough space
    }

    while (v_idx >= this->children_v.size()) {
      this->children_v.push_back(std::vector<pt::idx_t>()); // ensure children_v has enough space
    }

    this->vertices.push_back(v);
    v.set_idx(v_idx); // set the index of the vertex
    return v_idx;
  }

  void add_edge(pt::idx_t parent, pt::idx_t child) {
    while (child >= this->parent_v.size()) { this->parent_v.push_back(INVALID_ID); }
    this->parent_v[child] = parent;
    while (parent >= this->children_v.size()) { this->children_v.push_back(std::vector<pt::idx_t>()); }
    this->children_v[parent].push_back(child);
  }

  void del_edge(pt::idx_t parent, pt::idx_t child) {
    this->parent_v[child] = INVALID_IDX;

    std::vector<pt::idx_t>&children = this->children_v[parent];
    auto it = std::find(children.begin(), children.end(), child);
    if (it != children.end()) {
      children.erase(it);
    }
  }

  // ----
  // misc
  // ----
  void print_dot() const {
  std::cout << std::format(
    "graph G {{\n"
    "\trankdir = TD;\n"
    "\tnode[shape = circle];\n"
    "\tedge [arrowhead=vee];\n"
  );

  // print vertices
  for (pt::idx_t i{}; i < this->vtx_count(); i++) {
    std::cout << std::format("\t{} [label=\"{}\"];\n", i, this->get_vertex(i).as_str());
  }

  // print edges
  for (pt::idx_t i{}; i < this->vtx_count(); i++) {
    for (pt::idx_t c : this->get_children(i)) {
      std::cout << std::format("\t{} -- {};\n", i, c);
    }
  }

  std::cout << "}" << std::endl;
}

};


} // namespace povu::tree


#endif
