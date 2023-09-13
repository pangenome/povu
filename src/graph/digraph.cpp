#include <algorithm>
#include <cstddef>
#include <cstring>
#include <format>
#include <iostream>
#include <set>
#include <stack>
#include <string>
#include <vector>
#include <functional>
// #include <cstring>

#include "./digraph.hpp"

namespace hg = handlegraph;

namespace digraph {


// impl < for path_t
// bool operator<(const digraph::path_t& lhs, const digraph::path_t& rhs) {
//   return std::make_tuple(lhs.name, lhs.id) < std::make_tuple(rhs.name, rhs.id);
// }

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

// TODO: decouple internal graph rep with handlegraph ideas and names
/*
 * Vertex
 * ------
 */
Vertex::Vertex() :
  o(std::set<Edge>{}),
  i(std::set<Edge>{}),
  seq(std::string{}),
  handle(std::string{}),
  paths(std::set<std::size_t>{})
  {};

Vertex::Vertex(const std::string& sequence, const handlegraph::nid_t& id) :
  o(std::set<Edge>{}),
  i(std::set<Edge>{}),
  seq(sequence),
  handle(std::to_string(id)),
  paths(std::set<std::size_t>{})
{};

Vertex::Vertex(const std::string& sequence) :
  o(std::set<Edge>{}),
  i(std::set<Edge>{}),
  seq(sequence),
  handle(std::string{}),
  paths(std::set<std::size_t>{})
{};


std::set<Edge> const& Vertex::out() const { return this->o; }
std::set<Edge> const& Vertex::in() const { return this->i; }

std::set<Edge>* Vertex::out_mut() { return &this->o; }
std::set<Edge>* Vertex::in_mut() { return &this->i; }

//void Vertex::set_seq(std::string&& s) { this->seq = s;  }
//void Vertex::set_handle(std::string&& h) { this->handle = h; }

std::string const& Vertex::get_seq() const { return this->seq; }
std::string const& Vertex::get_handle() const { return this->handle; }

int Vertex::set_path(std::size_t p_id) {
  if (this->handle == "" || this->seq == "" ) { return 1; }

  this->paths.insert(p_id);
  return 0;
}

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

bool DiGraph::has_node(handlegraph::nid_t node_id) const {
		return node_id < this->size();
}

handlegraph::handle_t DiGraph::get_handle(const handlegraph::nid_t& node_id,
										  bool is_reverse) const {
  handlegraph::handle_t h;

  std::snprintf(h.data, sizeof(h.data), "%lld", node_id);

  //handle.data = std::to_string(node_id).c_str();

  return h;
}

handlegraph::nid_t DiGraph::get_id(const handlegraph::handle_t& handle) const {
  return std::stoi(handle.data);
}

bool DiGraph::get_is_reverse(const handlegraph::handle_t& handle) const {
  return false;
}

handlegraph::handle_t DiGraph::flip(const handlegraph::handle_t& handle) const {
	return handle;
}

size_t DiGraph::get_length(const handlegraph::handle_t& handle) const {
  return 0;
}

std::string DiGraph::get_sequence(const handlegraph::handle_t& handle) const {
	return "";
}

std::size_t DiGraph::get_node_count() const {
	return this->size();
}

handlegraph::nid_t DiGraph::min_node_id() const {
	return 0;
}

handlegraph::nid_t DiGraph::max_node_id() const {
	return 0;
}

bool DiGraph::follow_edges_impl(const handlegraph::handle_t& handle,
					   bool go_left,
					   const std::function<bool(const handlegraph::handle_t&)>& iteratee) const {
  return false;
}


bool DiGraph::for_each_handle_impl(const std::function<bool(const handlegraph::handle_t&)>& iteratee,
						  bool parallel) const {
		return false;
  }

// MutableHandleGraph
// ------------------
handlegraph::handle_t DiGraph::create_handle(const std::string& sequence) {
  handlegraph::handle_t h;

  //std::cout << "creating: " << sequence << this->size() << std::endl;
  std::snprintf(h.data, sizeof(h.data), "%ld", this->size());

  //std::cout << "creating: " << sequence << " size "<< this->size() << std::endl;
  this->adj.push_back(Vertex(sequence, this->size()));

  //std::cout << "creating: " << sequence << " size "<< this->size() << std::endl;

  //std::cout << "pushed\n";
  return h;
}

handlegraph::handle_t DiGraph::create_handle(const std::string& sequence, const handlegraph::nid_t& id) {
  std::size_t id_ = id;

  if (id < this->size()) {
	throw std::invalid_argument("id is less than size");
  }

  std::cout << "creating: " << this->size() << " " << id_ << std::endl;

  // pad with empty vertices until we reach the id
  for (std::size_t i = this->size(); i < id_; i++) {
	this->adj.push_back(Vertex());
  }

  return create_handle(sequence);
}

void DiGraph::create_edge(const handlegraph::handle_t& left, const handlegraph::handle_t& right) {
  this->add_edge(std::stoll(left.data), std::stoll(right.data));
}

handlegraph::path_handle_t DiGraph::create_path_handle(const std::string& name,
											  bool is_circular) {
  handlegraph::path_handle_t h;
  std::size_t path_id = this->paths.size();
  this->paths.push_back(path_t({name, path_id}));

  strncpy(h.data, std::to_string(path_id).c_str(), sizeof(h.data));

  //std::snprintf(h.data, sizeof(h.data), "%uld", this->size());
  return h;
}

handlegraph::path_handle_t
DiGraph::rename_path(const handlegraph::path_handle_t& path_handle,
					 const std::string& new_name) {
  // extract path id from path_handle
  std::size_t path_id = std::stoll(path_handle.data);
  this->paths[path_id].name = new_name;
  return path_handle;
}



//void DiGraph::add_path(const std::string& name) {
//  this->paths.insert(name);
//}

// constructor(s)
// --------------
DiGraph::DiGraph()
  : adj(std::vector<Vertex>{}),
	paths(std::vector<path_t>{})
{}

DiGraph::DiGraph(std::size_t size)
  : adj(std::vector<Vertex>{}),
	paths(std::vector<path_t>{})
{
  adj.reserve(size);
}

DiGraph::DiGraph(std::size_t size, std::size_t path_count)
  : adj(std::vector<Vertex>{}),
	paths(std::vector<path_t>{})
{
  adj.reserve(size);
  paths.reserve(path_count);
}

DiGraph::DiGraph(std::set<std::size_t>&& start_nodes, std::set<std::size_t>&& stop_nodes)
  : adj(std::vector<Vertex>{}),
	start_nodes(std::move(start_nodes)),
	end_nodes(std::move(stop_nodes)),
	paths(std::vector<path_t>{})
{}

// getters
// -------
void DiGraph::add_start_node(std::size_t idx) {
  this->start_nodes.insert(idx);
}
void DiGraph::add_stop_node(std::size_t idx) {
  // TODO: confirm not out of range
  this->end_nodes.insert(idx);
}

void DiGraph::compute_start_nodes() {
  for (std::size_t i = 0; i < this->size(); i++) {
	if (this->get_vertex(i).in().empty() && !this->get_vertex(i).out().empty()) {
	  this->add_start_node(i);
	}
  }
}

void DiGraph::compute_stop_nodes() {
  for (std::size_t i = 0; i < this->size(); i++) {
	if (this->get_vertex(i).out().empty() && !this->get_vertex(i).in().empty()) {
	  this->add_stop_node(i);
	}
  }
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
	// no throw an exception if the graph is too small
  for (std::size_t pos{size}; pos <= max; pos++) {
	// throw an out of range exception if the graph is too small
	throw std::out_of_range("graph is too small");
	//std::cout << "adding empty vertex\n";

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
