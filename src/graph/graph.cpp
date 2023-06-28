#include <algorithm>
#include <cstddef>
#include <format>
#include <iostream>
#include <limits>
#include <set>
#include <stack>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>

#include "./graph.hpp"



// undirected graph
namespace u_graph {

/*
 * Vertex
 */
Vertex::Vertex() : adj_vertices(std::set<std::size_t>{}) { }
void Vertex::add_adjacent_vertex(std::size_t vertex) {
  this->adj_vertices.insert(vertex);
}


std::set<std::size_t> const& Vertex::get_adjacent_vertices() const {
  return adj_vertices;
}

void Vertex::del_adjacent_vertex(std::size_t vertex) {
  this->adj_vertices.erase(vertex);
}


// at least 2
CFG::CFG(std::size_t initial_len) {
  std::vector<Vertex> d(initial_len, Vertex());
  this->adj_list = std::move(d);
  this->adj_list[0].add_adjacent_vertex(1);
  this->adj_list[1].add_adjacent_vertex(0);
}

std::size_t CFG::size_internal() {
  return this->adj_list.size();
}

std::size_t CFG::size() {
  return this->size_internal() - 2;
}

std::size_t CFG::stop_node_internal_idx() {
  return this->size_internal() - 1;
}

Vertex& CFG::start_node_internal() {
  return this->adj_list.at(this->start_node_id);
}

Vertex& CFG::stop_node_internal() {
  return this->adj_list.at(this->stop_node_internal_idx());
}


Vertex& CFG::get_vertex_mut(std::size_t vertex) {
  return this->adj_list.at(++vertex);
}

Vertex const& CFG::get_vertex(std::size_t vertex) const {
  return this->adj_list.at(++vertex);
}

Vertex const& CFG::get_vertex_internal(std::size_t vertex) const {
  return this->adj_list.at(vertex);
}

// will break when called with internal vertex indexes
std::set<std::size_t>const& CFG::get_adjacent_vertices(std::size_t vertex) const {
  return this->get_vertex(vertex).get_adjacent_vertices();
}

// will break when called with internal vertex indexes
std::set<std::size_t> const& CFG::get_adjacent_vertices_internal(std::size_t vertex) const {
  return this->get_vertex_internal(vertex).get_adjacent_vertices();
}

void CFG::set_start_node(std::size_t vertex) {
  this->start_node_internal().add_adjacent_vertex(vertex+1);
  this->get_vertex_mut(vertex).add_adjacent_vertex(this->start_node_id);
}

void CFG::set_stop_node(std::size_t vertex) {
  this->get_vertex_mut(vertex).add_adjacent_vertex(this->stop_node_internal_idx());
  this->stop_node_internal().add_adjacent_vertex(++vertex);
}

void CFG::add_edge(std::size_t n1, std::size_t n2) {
  ++n1; ++n2;
  // move the last node back if n2 or n1 reaches it
  // TODO: what about self loops i.e n1 == n2?
  std::size_t size = this->size();
  // TODO: update current start and stop nodes
  if (size < std::max(n1, n2)) {
    // delete edge from start to final
    // is out of range
    this->start_node_internal().del_adjacent_vertex(this->size_internal() - 1);

    std::size_t new_internal_size = this->size_internal() + std::max(n1, n2) - size;
    this->adj_list.resize(new_internal_size, Vertex());

    std::size_t stop_node_idx = size + 1;
    std::swap(this->adj_list[stop_node_idx], this->adj_list[new_internal_size - 1]);
    // add edge from start to new final node
    this->start_node_internal().add_adjacent_vertex(new_internal_size - 1);
  }

  this->adj_list[n1].add_adjacent_vertex(n2);
  this->adj_list[n2].add_adjacent_vertex(n1);
}

spanning_tree::Tree CFG::compute_spanning_tree() {
  spanning_tree::Tree t = spanning_tree::Tree(this->size_internal());

  //std::set<std::size_t> explored;
  std::set<std::size_t> seen;
  std::stack<std::size_t> visited;

  std::size_t current_vertex{this->start_node_id};
  visited.push(current_vertex);

  std::size_t counter{0};

  while (!visited.empty()) {
    current_vertex = visited.top();

    if (!seen.count(current_vertex)) {
      t.set_dfs_num(current_vertex, counter);
      t.set_sort(counter, current_vertex);
      ++counter;
    }

    seen.insert(current_vertex);

    // TODO: simplify below for loop
    // - replace f with not_explored
    // bool not_explored{false}; // the current vertex has not been explored
    bool f{false};
    Vertex const& v =  this->get_vertex_internal(current_vertex);

    // TODO: better condition here for speedup
    for (auto a : v.get_adjacent_vertices()) {
      if (seen.find(a) == seen.end()) {
        t.add_tree_edge(current_vertex, a);
        visited.push(a);
        f = true;
        break;
      }
      else if (
        !t.is_root(current_vertex) &&
        t.get_parent(current_vertex) != a &&
        !t.has_child(current_vertex, a) &&
        !t.has_ibe(current_vertex, a) &&
        !t.has_obe(current_vertex, a)
      ) {
        // TODO: why the has child and not parent test?
         t.add_be(current_vertex, a);
      }
    }


    if (!f) {
      visited.pop();
    }
  }

  return t;
}

  
void CFG::print_dot() {
  std::cout << std::format(
    "graph G {{\n"
    "\trankdir = TB;\n"
    "\tnode[shape = circle];\n"
    "\tedge [arrowhead=vee];\n"
    "\t{} [label=\"S\", color=\"green\"];\n"
    "\t{} [label=\"E\", color=\"green\"];\n",
    this->start_node_id,
    this->stop_node_internal_idx());

  std::set<Edge> reported{};


  for (std::size_t i{}; i < this->size_internal(); i++) {
    for (auto o : this->get_adjacent_vertices_internal(i)) {
      if (reported.find(Edge(i, o)) != reported.end()) { continue; }

      reported.insert(Edge(i,o));
      std::cout << std::format("\t{} -- {};\n", i, o);
    }
  }

  std::cout << "}" << std::endl;
}

}; // namespace u_graph
