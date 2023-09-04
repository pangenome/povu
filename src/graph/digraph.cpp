#include <algorithm>
#include <cstddef>
#include <format>
#include <iostream>
#include <set>
#include <stack>
#include <vector>


#include "./digraph.hpp"

namespace digraph {
/*
 * Edge
 * ----
 */
Edge::Edge() : frm(0), t(0), c(core::color::black) {};
Edge::Edge(std::size_t frm, std::size_t to, core::color c) : frm(frm), t(to), c(c) {};

// implement operator< for Edge
bool operator<(const Edge& lhs, const Edge& rhs) {
  return std::make_tuple(lhs.from(), lhs.to(), lhs.get_color()) <
    std::make_tuple(rhs.from(), rhs.to(), rhs.get_color());
}

std::size_t Edge::to() const { return this->t; }
std::size_t Edge::from() const { return this->frm; }
core::color Edge::get_color() const { return this->c; }
bool Edge::is_black() const { return this->c == core::color::black; }

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
void Vertex::add_out(std::size_t self_idx, std::size_t to_idx, core::color c) {
  this->o.insert(Edge(self_idx, to_idx, c));
};
void Vertex::add_in(std::size_t from_idx, std::size_t self_idx, core::color c) {
  this->i.insert(Edge(from_idx, self_idx, c));
};
bool Vertex::is_leaf() const { return this->out().empty(); }

/*
 * DiGraph
 * -------
 */

// handlegraph

handlegraph::handle_t DiGraph::create_handle(const std::string& sequence) {
  this->create_handle(sequence, this->adj.size());
}

handlegraph::handle_t DiGraph::create_handle(const std::string& sequence, const handlegraph::nid_t& id) {
  // we expect id to be equal to the end of the vector

  handlegraph::handle_t h;
  
  if (id == this->size()) {
	this->adj.push_back(Vertex());
  }
  else if (id > this->adj.size()) {
	//this->adj.resize(id + 1);
  }
  else {
	// TODO: throw an error because this edge exists
	// TODO: check if already exists
  }
  //this->adj.push_back(Vertex());
  //Vertex()
  //this->create_handle(sequence, this->adj.size());

  return h;
}

void DiGraph::create_edge(const handlegraph::handle_t& left, const handlegraph::handle_t& right) {
  std::size_t l_value = std::stoul(left.data);
  std::size_t r_value = std::stoul(right.data);

  this->add_edge(l_value, r_value, core::color::black);
}

// FIXME: does this work in our context?
handlegraph::handle_t DiGraph::apply_orientation(const handlegraph::handle_t& handle){
  handlegraph::handle_t h;
  return h;
}

std::vector<handlegraph::handle_t>
DiGraph::divide_handle(const handlegraph::handle_t& handle, const std::vector<std::size_t>& offsets) {
	std::vector<handlegraph::handle_t> handles;
	return handles;
}

void DiGraph::optimize(bool allow_id_reassignment) {
	return;
}

bool DiGraph::apply_ordering(const std::vector<handlegraph::handle_t>& order, bool compact_ids) {
  	return false;
}

void DiGraph::set_id_increment(const handlegraph::nid_t& min_id) {
	return;
}

void DiGraph::increment_node_ids(handlegraph::nid_t increment) {
	
}
  
void DiGraph::increment_node_ids(long increment) {}
    
void DiGraph::reassign_node_ids(const std::function<handlegraph::nid_t(const handlegraph::nid_t&)>& get_new_id) {}
  
// constructor(s)
// --------------
DiGraph::DiGraph() : adj(std::vector<Vertex>{}) {};
DiGraph::DiGraph(std::size_t size) : adj(std::vector<Vertex>{}) {
  adj.reserve(size);
}
DiGraph::DiGraph(std::set<std::size_t>&& start_nodes, std::set<std::size_t>&& stop_nodes)
  : adj(std::vector<Vertex>{}),
    start_nodes(std::move(start_nodes)),
    end_nodes(std::move(stop_nodes)) {};

// getters
// -------
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

void DiGraph::add_edge(std::size_t from, std::size_t to, core::color c) {

  std::size_t size = this->size();
  std::size_t max = std::max(from, to);

  // if the graph is not big enough, add empty vertices
  for (std::size_t pos{size}; pos <= max; pos++) {
    this->adj.push_back(Vertex());
  }

  this->adj[from].add_out(from, to, c);
  this->adj[to].add_in(from, to, c);
}

/**
 * If an edge has more than 1 in coming AND more than
 * 1 outgoing, then split it into two nodes.
 *
 * To bi-edge a graph, we need to:
 * 1. Find all nodes with more than 1 incoming or outgoing edge
 * 2. For each node, split it into two nodes and add a gray edge between them
 * incoming nodes are added to the first node, outgoing to the second
 * this is because the gray edge is directed from n to n+1
 * a node n is split into two nodes n and n+1 where n is the first (or left)
 * node & n+1 is the second (or right) node
 */
void DiGraph::biedge() {
 
  for (std::size_t idx{}; idx < this->adj.size(); idx++) {
    Vertex const& v = this->get_vertex(idx);

    if (
      (v.in().size() > 1 && v.out().size() > 1) ||
      ((v.out().size() > 1) && v.in().begin()->is_black()) ||
      ((v.in().size() > 1) && v.out().begin()->is_black())
      )
    {

      // duplicate the node (v_0) at idx
      Vertex v_0 = this->adj[idx];
      this->adj.insert(this->adj.begin() + idx, v_0);

      Vertex& v1 = this->get_vertex_mut(idx);
      Vertex& v2 = this->get_vertex_mut(idx+1);

      v1.out_mut()->clear();
      v2.in_mut()->clear();

      this->add_edge(idx, idx+1, core::color::gray);

      // increment all affected edges
      for (std::size_t j{}; j < this->adj.size(); j++) {

        if (j == idx) {  continue; }

        Vertex* v = &this->adj[j];
        std::set<Edge>* in = v->in_mut();
        std::set<Edge>* out = v->out_mut();

        for (auto& e : *out) {
          if (e.to() > idx) { const_cast<Edge&>(e).set_to(e.to() + 1); }
        }

        for (auto& e : *in) {
          if (e.from() > idx) { const_cast<Edge&>(e).set_from(e.from() + 1); }
        }
      }

      // increment stop nodes
      for (auto& n : this->end_nodes) {
        if (n > idx) {
          const_cast<std::size_t&>(n) = n + 1;
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
      std::cout << std::format("\t{} -> {} [color=\"{}\"];\n",
                               i,
                               o.to(),
                               (o.get_color() == core::color::black ? "black" : "gray")
        );
    }
  }
  std::cout << "}" << std::endl;
}
} // namespace digraph
