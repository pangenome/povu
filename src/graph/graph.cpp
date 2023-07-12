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
#include "./digraph.hpp"


// undirected graph
namespace u_graph {

/*
 * Vertex
 * ======
 */
  
// constructor(s)
// --------------
Vertex::Vertex() : adj_vertices(std::set<std::size_t>{}) { }


// getters
// -------
std::set<std::size_t> const& Vertex::get_adjacent_vertices() const {
  return adj_vertices;
}
  
// setters
// -------
void Vertex::add_adjacent_vertex(std::size_t vertex) {
  this->adj_vertices.insert(vertex);
}

void Vertex::del_adjacent_vertex(std::size_t vertex) {
  this->adj_vertices.erase(vertex);
}

/*
 * Flow Graph
 * ==========
 */

// at least 2
CFG::CFG(std::size_t initial_len) {
  std::vector<Vertex> d(initial_len, Vertex());
  this->adj_list = std::move(d);
  std::size_t edge_idx = this->edges.size();
  
  this->edges.push_back(u_graph::Edge{0, initial_len - 1});

  this->adj_list[0].add_adjacent_vertex(edge_idx);
  this->adj_list[initial_len - 1].add_adjacent_vertex(edge_idx);
}
  
CFG::CFG(digraph::DiGraph const& di_graph) {
  std::size_t size = di_graph.size();

  // initialize the flow graph
  std::vector<Vertex> d(size+2, Vertex());
  this->adj_list = std::move(d);


  std::size_t edge_idx = this->edges.size();
  
  this->edges.push_back(u_graph::Edge{0, size + 1});

  // connect start to end
  this->adj_list[0].add_adjacent_vertex(edge_idx);
  this->adj_list[size + 1].add_adjacent_vertex(edge_idx);

  // connect all digraph start nodes to flow graph start node
  for (auto const& start_node : di_graph.starts()) {
    edge_idx = this->edges.size();
    this->edges.push_back(u_graph::Edge{0, start_node+1});
    this->adj_list[0].add_adjacent_vertex(edge_idx);
    // no inc because zero index
    this->adj_list[start_node].add_adjacent_vertex(edge_idx);
  }

  for (std::size_t i{}; i < size; ++i) {
    for (auto const& edge : di_graph.get_vertex(i).out()) {
      this->add_edge(i, edge.to(), edge.get_color());
    }
  }


  // connect all digraph end nodes to stop node
  for (auto const& stop_node : di_graph.stops()) {
    edge_idx = this->edges.size();
    
    this->edges.push_back(u_graph::Edge{stop_node+1, this->stop_node_internal_idx()});

    this->adj_list[this->stop_node_internal_idx()].add_adjacent_vertex(edge_idx);
    this->adj_list[stop_node+1].add_adjacent_vertex(edge_idx);
    //this->set_stop_node(stop_node);
  }
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

void CFG::add_edge(std::size_t n1, std::size_t n2, core::color c) {
  ++n1; ++n2;
  // move the last node back if n2 or n1 reaches it
  // TODO: what about self loops i.e n1 == n2?
  std::size_t size = this->size();
  std::size_t greater = std::max(n1, n2);

  // TODO: update current start and stop nodes
  if (size < greater) {
    std::size_t new_internal_size = this->size_internal() + greater - size;
    
    // delete edge from start to final
    // is out of range
    for (auto edge_idx: this->start_node_internal().get_adjacent_vertices()) {

      Edge& edge = this->edges.at(edge_idx);
      
      // TODO: check if stop_node is left
      if (edge.right() == this->stop_node_internal_idx()) {
        edge.set_right(new_internal_size - 1);
        //this->start_node_internal().del_adjacent_vertex(e);
        break;
      }
    }
    
    //this->start_node_internal().del_adjacent_vertex(this->size_internal() - 1);

    // std::size_t new_internal_size = this->size_internal() + greater - size;
    
    this->adj_list.resize(new_internal_size, Vertex());

    std::size_t curr_stop_node_idx = size + 1;
    std::swap(this->adj_list[curr_stop_node_idx],
              this->adj_list[new_internal_size - 1]);
    
    // add edge from start to new final node
    // this->start_node_internal().add_adjacent_vertex(new_internal_size - 1);
  }

  std::size_t edge_idx = this->edges.size();  
  this->edges.push_back(u_graph::Edge{n1, n2, c});
  this->adj_list[n1].add_adjacent_vertex(edge_idx);
  this->adj_list[n2].add_adjacent_vertex(edge_idx);
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
    for (auto edge_idx : v.get_adjacent_vertices()) {
      Edge& edge = this->edges.at(edge_idx);
      std::size_t a = edge.left() == current_vertex ? edge.right() : edge.left();
      
      if (seen.find(a) == seen.end()) {
        t.add_tree_edge(current_vertex, a, edge.get_color());
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
        //std::cout << "adding back edge: " << current_vertex << " -> " << a << std::endl;
        t.add_be(current_vertex, a, false, edge.get_color());
      }
    }

    if (!f) { visited.pop(); }
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

  //std::set<Edge> reported{};
  std::set<std::size_t> reported{};

  for (std::size_t i{}; i < this->size_internal(); i++) {
    for (std::size_t  edge_idx : this->get_adjacent_vertices_internal(i)) {
      if (reported.count(edge_idx)) { continue; }

      reported.insert(edge_idx);

      Edge& edge = this->edges.at(edge_idx);
      
      std::size_t to = edge.right() == i ? edge.left() : edge.right();
      
      std::cout << std::format("\t{} -- {}  [color=\"{}\"];\n",
                               i, to,
                               (this->edges.at(edge_idx).get_color() == core::color::black ? "black" : "gray"));
    }
  }

  std::cout << "}" << std::endl;
}

}; // namespace u_graph
