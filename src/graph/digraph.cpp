#include <algorithm>
#include <cstddef>
#include <format>
#include <set>
#include <stack>
#include <vector>


#include "./digraph.hpp"

namespace digraph {
/*
 * Edge
 * ----
 */
Edge::Edge() : frm(0), t(0), c(colour::black) {};
Edge::Edge(std::size_t frm, std::size_t to, colour c) : frm(frm), t(to), c(c) {};

// implement operator< for Edge
bool operator<(const Edge& lhs, const Edge& rhs) {
  return std::make_tuple(lhs.from(), lhs.to(), lhs.get_colour()) <
    std::make_tuple(rhs.from(), rhs.to(), rhs.get_colour());
}

std::size_t Edge::to() const { return this->t; }  
std::size_t Edge::from() const { return this->frm; }
colour Edge::get_colour() const { return this->c; }

void Edge::set_to(std::size_t t) { this->t = t; }
void Edge::set_from(std::size_t f) { this->frm = f; }  

  
/*
 * Vertex
 * ------
 */
Vertex::Vertex() : o(std::set<Edge>{}), i(std::set<Edge>{}) {};

std::set<Edge> const& Vertex::out() const { return this->o; }
std::set<Edge> const& Vertex::in() const { return this->i; }

std::set<Edge>* Vertex::out_mut() { return &this->o; }
std::set<Edge>* Vertex::in_mut() { return &this->i; }
  
  
// TODO: not use zero
void Vertex::add_out(std::size_t self_idx, std::size_t to_idx) {
  this->o.insert(Edge(self_idx, to_idx));
};
void Vertex::add_in(std::size_t from_idx, std::size_t self_idx) {
  this->i.insert(Edge(from_idx, self_idx));
};
bool Vertex::is_leaf() const { return this->out().empty(); }

/*
 * DiGraph
 * -------
 */

DiGraph::DiGraph() : adj(std::vector<Vertex>{}) {};
DiGraph::DiGraph(std::size_t size) : adj(std::vector<Vertex>{}) {
  adj.reserve(size);
} 
DiGraph::DiGraph(std::set<std::size_t>&& start_nodes, std::set<std::size_t>&& stop_nodes)
  : adj(std::vector<Vertex>{}),
    start_nodes(std::move(start_nodes)),
    end_nodes(std::move(stop_nodes)) {};

void DiGraph::add_start_node(std::size_t idx) {
  this->start_nodes.insert(idx);
}
void DiGraph::add_stop_node(std::size_t idx) {
  // TODO: confirm not out of range
  this->end_nodes.insert(idx);
}

Vertex const& DiGraph::get_vertex(std::size_t idx) const {
  return this->adj.at(idx);
}

Vertex& DiGraph::get_vertex_mut(std::size_t idx) {
  return this->adj.at(idx);
}

std::set<std::size_t> const& DiGraph::starts() const {
  return this->start_nodes;
}

std::set<std::size_t> const& DiGraph::stops() const {
  return this->end_nodes;
}

std::size_t DiGraph::size() const { return this->adj.size(); }

void DiGraph::add_edge(std::size_t from, std::size_t to) {
  
  std::size_t size = this->size();
  std::size_t max = std::max(from, to);

  // if the graph is not big enough, add empty vertices
  for (std::size_t pos{size}; pos <= max; pos++) {
    this->adj.push_back(Vertex());
  }

  this->adj[from].add_out(from, to);
  this->adj[to].add_in(from, to);
}

/**
 * If an edge has more than 1 in coming AND more than
 * 1 outgoing, then split it into two nodes.
 */
void DiGraph::biedge() {
  for (std::size_t idx{}; idx < this->adj.size(); idx++) {
    Vertex const& v = this->get_vertex(idx);
    
    if (v.out().size() > 1 && v.in().size() > 1) {
      
      Vertex v_0 = this->adj[idx];
      this->adj.insert(this->adj.begin() + idx, v_0);

      Vertex& v1 = this->get_vertex_mut(idx);
      Vertex& v2 = this->get_vertex_mut(idx+1);

      v1.out_mut()->clear();
      v2.in_mut()->clear();
      
      this->add_edge(idx, idx+1);

      // increment all edges after idx
      for (std::size_t j{}; j < this->adj.size(); j++) {

        if (j == idx) {  continue; }

        Vertex* v = &this->adj[j];
        std::set<Edge>* out = v->out_mut();
        std::set<Edge>* in = v->in_mut();

        for (auto& e : *out) {
          if (e.to() > idx) {
            const_cast<Edge&>(e).set_to(e.to() + 1);
          }
        }

        for (auto& e : *in) {
          if (e.from() > idx) {
            const_cast<Edge&>(e).set_from(e.from() + 1);
          }
        }
      }
    }
  }
}
  
void DiGraph::print_dot() {
  std::cout << std::format(
    "digraph G {{\n"
    "\trankdir = TB;\n"
    "\tnode[shape = circle];\n"
    "\tedge [arrowhead=vee];\n"
  );

  for (auto n: this->starts()) {
    std::cout << std::format("\t{} [color=\"green\"];\n", n);
  }

  for (auto n: this->stops()) {
    std::cout << std::format("\t{} [color=\"blue\"];\n", n);
  }

  for (std::size_t i{}; i < this->size(); i++) {
    // for each outgoing vertex
    for (auto o : this->get_vertex(i).out()) {
      std::cout << std::format("\t{} -> {};\n", i, o.to());
    }
  }
  std::cout << "}" << std::endl;
}
} // namespace digraph
