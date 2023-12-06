#include <cstddef>
#include <format>
#include <iostream>
#include <string>
#include <vector>
#include <stack>

#include "./biedged.hpp"
#include "./bidirected.hpp"

namespace biedged {

/* VertexType
 */

// implement << operator for VertexType
std::ostream& operator<<(std::ostream& os, const VertexType& vt) {
	switch (vt) {
	case VertexType::l:
	os << "+";
	break;
	case VertexType::r:
	os << "-";
	break;
	default:
	os << "*";
	break;
	}

	return os;
}

/*
  Edge
*/


Edge::Edge(std::size_t v1, VertexType v1_type,
		   std::size_t v2, VertexType v2_type,
		   core::color c)
  : v1_idx(v1), v1_type(v1_type), v2_idx(v2), v2_type(v2_type), c(c)
{

  if (v1 == v2) {
	throw std::invalid_argument(
	  std::format("[povu::graph::biedged::Edge] Self-loops are not allowed {} {}", v1, v2));
  }

  if (c == core::color::black) {
	// throw an argument error if the edge is black and no label is provided
	throw std::invalid_argument("[povu::graph::biedged::Edge] Black edges must have a label");
  }

  this->label = std::string();

}

Edge::Edge(std::size_t v1, VertexType v1_type,
		   std::size_t v2, VertexType v2_type,
		   core::color c, std::string label)
  :  v1_idx(v1), v1_type(v1_type), v2_idx(v2), v2_type(v2_type), c(c), label(label)
{

  if (v1 == v2) {
	throw std::invalid_argument("[povu::graph::biedged::Edge] Self-loops are not allowed");
  }

  if (c == core::color::gray) {
	// throw an argument error if the edge is gray and a label is provided
	throw std::invalid_argument("[povu::graph::biedged::Edge] Gray edges cannot have a label");
  }
}


std::size_t Edge::get_v1_idx() const {
	return this->v1_idx;
}

std::size_t Edge::get_v2_idx() const {
	return this->v2_idx;
}

core::color Edge::get_color() const {
  return this->c;
}

const std::string& Edge::get_label() const {
	return this->label;
}

void Edge::set_v1_idx(std::size_t i) {
	this->v1_idx = i;
}

void Edge::set_v2_idx(std::size_t i) {
	this->v2_idx = i;
}

std::ostream& operator<<(std::ostream& os, const Edge& e) {
	os << "Edge(" << e.v1_idx << ", "
	   << (e.v1_type == VertexType::l ? "-" : "+")
	   << e.v2_idx << ", "
	   << (e.v2_type == VertexType::l ? "-" : "+")
	   << e.c << ")";
	return os;
}

/*
  Vertex
*/
Vertex::Vertex() :
  black_edge(0),
  grey_edges(std::set<std::size_t>()),
  paths(std::vector<PathInfo>()),
  type(VertexType::l),
  handle(std::string()),
  is_reversed_(false),
  vertex_idx(0)
{}

Vertex::Vertex(const std::string& id, std::size_t vertex_idx, VertexType vertex_type) :
	black_edge(0),
	grey_edges(std::set<std::size_t>()),
	paths(std::vector<PathInfo>()),
	type(vertex_type),
	handle(id),
	is_reversed_(false),
	vertex_idx(vertex_idx)
{}


const std::string& Vertex::get_handle() const {
	return this->handle;
}

std::set<std::size_t> Vertex::get_grey_edges() const {
	return this->grey_edges;
}

std::size_t Vertex::get_black_edge() const {
	return this->black_edge;
}

std::size_t Vertex::get_vertex_idx() const {
  return this->vertex_idx;
}

void Vertex::add_edge(std::size_t edge_idx, core::color c) {
  if (c == core::color::black) {
	black_edge = edge_idx;
  } else {
	grey_edges.insert(edge_idx);
  }
}

VertexType Vertex::get_type() const {
	return type;
}

void Vertex::set_vertex_idx(std::size_t i) {
	this->vertex_idx = i;
}

/*
  BVariationGraph
 */

void BVariationGraph::add_edge(std::size_t v1,VertexType v1_type,
							   std::size_t v2, VertexType v2_type,
							   core::color c) {
  edges.push_back(Edge(v1, v1_type, v2, v2_type, c));
}

void BVariationGraph::add_edge(std::size_t v1, VertexType v1_type,
							   std::size_t v2, VertexType v2_type,
							   core::color c, std::string label) {
  edges.push_back(Edge(v1, v1_type, v2, v2_type, c, label));
}

void BVariationGraph::add_vertex(std::string handle_str, VertexType vertex_type) {
  vertices.push_back(Vertex(handle_str, vertices.size(), vertex_type));
}

bool BVariationGraph::replace_vertex(std::size_t vertex_idx, std::string handle_str, VertexType vertex_type) {
	if (this->get_vertex(vertex_idx).get_handle() != "") {
		return false;
	}

	this->get_vertex_mut(vertex_idx) = Vertex(handle_str, vertex_idx, vertex_type);
	return true;
}

const biedged::Vertex& BVariationGraph::get_vertex(std::size_t i) const {
	return vertices.at(i);
}

std::size_t BVariationGraph::size() const {
	return vertices.size();
}

Vertex& BVariationGraph::get_vertex_mut(std::size_t i) {
	return vertices.at(i);
}

BVariationGraph::BVariationGraph(const bidirected::VariationGraph& g) {

  this->vertices =
	std::vector<biedged::Vertex>(g.size()*2, biedged::Vertex());

  auto new_l = [](std::size_t i) { return (2*(i+1)) - 2; };
  auto new_r = [](std::size_t i) { return (2*(i+1)) - 1; };

  /*

   */
  auto duplicate_vertices =
	[this, g](bidirected::Vertex const& current_vg_vertex, std::size_t i_l, std::size_t i_r) {

	  // --------------------
	  // add vertices + and -
	  // --------------------

	  {
		this->replace_vertex(i_l, current_vg_vertex.get_handle(), biedged::VertexType::l);
		this->replace_vertex(i_r, current_vg_vertex.get_handle(), biedged::VertexType::r);

	  }

	  // ----------------
	  // add a black edge
	  // ----------------

	  {
		this->add_edge(
		  i_l, VertexType::l,
		  i_r, VertexType::r,
		  core::color::black,
		  g.get_vertex(0).get_label());

		this->get_vertex_mut(i_l).add_edge(
		  this->edges.size() - 1, core::color::black);
		this->get_vertex_mut(i_r).add_edge(
		  this->edges.size() - 1, core::color::black);
	  }
	};

  std::set<unordered_pair> added_edges;

  /*
	add gray edges
  */
  auto do_gray_edges =
	[&](std::size_t v_idx, bidirected::Vertex const& current_vg_vertex, std::size_t i_l, std::size_t i_r) {


	  // add gray edges incident with the 5' (left/+) vertex
	  for (auto e_idx : current_vg_vertex.get_edges_l()) {

		bidirected::Edge e = g.get_edge(e_idx);

		std::size_t new_v1{};
		VertexType v1_type{};


		if (e.get_v1_idx() == v_idx){
		  new_v1 = e.get_v2_end() == bidirected::VertexEnd::l ?
			new_l(e.get_v2_idx()) :
			new_r(e.get_v2_idx());

		  v1_type =
			e.get_v2_end() == bidirected::VertexEnd::l ?
			VertexType::l : VertexType::r;
		} else {
		  new_v1 =
			e.get_v1_end() == bidirected::VertexEnd::l ?
			new_l(e.get_v1_idx()) :
			new_r(e.get_v1_idx());

		  v1_type =
			e.get_v1_end() == bidirected::VertexEnd::l ?
			VertexType::l : VertexType::r;
		}


		if (added_edges.count(unordered_pair(new_v1, i_l))) { continue; }

		added_edges.insert(unordered_pair(new_v1, i_l));

		this->add_edge(new_v1, v1_type,
					   i_l, VertexType::l,
					   core::color::gray);

		this->get_vertex_mut(new_v1).add_edge(
		  this->edges.size() - 1, core::color::gray);
		this->get_vertex_mut(i_l).add_edge(
		  this->edges.size() - 1, core::color::gray);
	  }

	  // add gray edges incident with the 3' (right/-) vertex
	  for (auto e_idx : current_vg_vertex.get_edges_r()) {

		bidirected::Edge e = g.get_edge(e_idx);

		std::size_t new_v2{};
		VertexType v2_type;

		if (e.get_v2_idx() == v_idx){
		  new_v2 =
			e.get_v1_end() == bidirected::VertexEnd::l ?
			new_l(e.get_v1_idx()) :
			new_r(e.get_v1_idx());

		  v2_type =
			e.get_v1_end() == bidirected::VertexEnd::l ?
			VertexType::l : VertexType::r;
		}
		else {
		  new_v2 =
			e.get_v2_end() == bidirected::VertexEnd::l ?
			new_l(e.get_v2_idx()) :
			new_r(e.get_v2_idx());

		  v2_type =
			e.get_v2_end() == bidirected::VertexEnd::l ? VertexType::l : VertexType::r;
		}

		if (added_edges.count(unordered_pair(i_r, new_v2))) { continue; }

		added_edges.insert(unordered_pair(i_r, new_v2));

		this->add_edge(i_r, VertexType::r,
					   new_v2, v2_type,
					   core::color::gray);

		this->get_vertex_mut(i_r).add_edge(
		  this->edges.size() - 1, core::color::gray);
		this->get_vertex_mut(new_v2).add_edge(
		  this->edges.size() - 1, core::color::gray);
	  }

	return;
  };


  for (std::size_t i = 0; i < g.size(); ++i) {
	duplicate_vertices(g.get_vertex(i), (2*(i+1)) - 2 , (2*(i+1)) - 1);
  }

  for (std::size_t i = 0; i < g.size(); ++i) {
	do_gray_edges(i, g.get_vertex(i),  (2*(i+1)) - 2 , (2*(i+1)) - 1);
  }

  // set start and stop vertices
  for (auto i : g.get_start_nodes()) {
	this->start_nodes.insert((2*(i+1)) - 2);
  }


  for (auto i : g.get_end_nodes()) {
	this->end_nodes.insert((2*(i+1)) - 1);
  }

}

void BVariationGraph::print_dot() const {
  //
  std::cout << "graph G {\n" <<
	"\trankdir=LR;\n" <<
	"\tnode [shape=circle];\n";

  //
  for (std::size_t v_idx{}; v_idx < this->vertices.size(); ++v_idx) {
	std::cout <<
	  std::format("\t{} [label=\"{}{}\"]\n",
				  v_idx,
				  this->get_vertex(v_idx).get_handle(),
				  (this->get_vertex(v_idx).get_type() == VertexType::l ? "+" : "-") );
  }

  //
  for (std::size_t e_idx{}; e_idx < this->edges.size() ; ++e_idx) {
	auto e = this->edges.at(e_idx);
	if (e.get_color() == core::color::black) {
	  std::cout <<
		std::format("\t{} -- {} [color=\"black\"; label=\"{}\"];\n",
					e.get_v1_idx(), e.get_v2_idx(), e.get_label());
	}
	else {
	  std::cout <<
		std::format("\t{} -- {} [color=\"gray\"];\n", e.get_v1_idx(), e.get_v2_idx());
	}
  }

  std::cout << "}" << std::endl;
}

void BVariationGraph::componetize() {
  this->vertices.insert(this->vertices.begin(), Vertex("d_s", 0, VertexType::dummy));
  this->vertices.push_back(Vertex("d_e", this->size(), VertexType::dummy));


  // ------------------------------
  // Update vertex indices in edges
  // ------------------------------
  for (std::size_t e_idx{}; e_idx < this->edges.size(); ++e_idx) {
	Edge& e = this->edges.at(e_idx);

	e.set_v1_idx(e.get_v1_idx() + 1);
	e.set_v2_idx(e.get_v2_idx() + 1);
  }


  // ---------------------------------
  // Update vertex indices in vertices
  // ---------------------------------
  for (std::size_t v_idx{1}; v_idx < this->vertices.size()-1; ++v_idx) {
	Vertex& v = this->get_vertex_mut(v_idx);
	v.set_vertex_idx(v.get_vertex_idx() + 1);
  }


  // -------------------------------------
  // connect dummy start to start vertices
  // -------------------------------------
  for (std::size_t s : this->start_nodes) {
	this->add_edge(
	  0, VertexType::dummy,
	  s+1, VertexType::l,
	  core::color::gray);

	this->get_vertex_mut(0)
	  .add_edge(this->edges.size() - 1, core::color::gray);
	this->get_vertex_mut(s+1)
	  .add_edge(this->edges.size() - 1, core::color::gray);
  }


  // -----------------------------------
  // connect dummy stop to stop vertices
  // -----------------------------------
  for (std::size_t s : this->end_nodes) {
	this->add_edge(
	  this->size() - 1, VertexType::dummy,
	  s+1, VertexType::l,
	  core::color::gray);

	this->get_vertex_mut(this->size() - 1)
	  .add_edge(this->edges.size() - 1, core::color::gray);
	this->get_vertex_mut(s+1)
	  .add_edge(this->edges.size() - 1, core::color::gray);
  }


  // ---------------------------------
  // connect dummy start to dummy stop
  // ---------------------------------
  this->add_edge(
	0, VertexType::dummy,
	this->size() - 1, VertexType::dummy,
	core::color::gray);

  this->get_vertex_mut(0)
	.add_edge(this->edges.size() - 1, core::color::gray);
  this->get_vertex_mut(this->size() -1)
	.add_edge(this->edges.size() - 1, core::color::gray);
}

spanning_tree::Tree BVariationGraph::compute_spanning_tree() const {

  spanning_tree::Tree t = spanning_tree::Tree(this->size());

  std::set<std::size_t> seen;
  std::stack<std::size_t> visited;

  std::size_t start_node_id{};

  std::size_t current_vertex{start_node_id};
  visited.push(current_vertex);

  std::size_t counter{0};

  while (!visited.empty()) {
	current_vertex = visited.top();

	//std::cout << "current_vertex: " << current_vertex << std::endl;
	
	if (!seen.count(current_vertex)) {
	  t.set_dfs_num(current_vertex, counter);
	  t.set_sort(counter, current_vertex);
	  t.set_sort_g(current_vertex, counter);
	  ++counter;
	}

	seen.insert(current_vertex);

	// TODO: simplify below for loop
	// - replace f with not_explored
	// bool not_explored{false}; // the current vertex has not been explored
	bool f{false};
	Vertex const& v =  this->get_vertex(current_vertex);

	std::set<size_t> adj_edges = v.get_grey_edges();
	if (v.get_type() != VertexType::dummy) { adj_edges.insert(v.get_black_edge()); }

	// target vertex, edge index
	std::vector<std::pair<std::size_t, std::size_t>> adj_vertices;
	for (std::size_t e_idx : adj_edges) {
	  Edge e = this->edges.at(e_idx);
	  std::size_t adj_v_idx = e.get_v1_idx() == current_vertex ? e.get_v2_idx() : e.get_v1_idx();
	  adj_vertices.push_back( std::make_pair(adj_v_idx, e_idx) );
	}

	std::sort(
	  adj_vertices.begin(),
	  adj_vertices.end(),
	  [](std::pair<std::size_t, std::size_t> const& a, std::pair<std::size_t, std::size_t> const& b) {
		return a.first < b.first;
	  });

	// TODO: better condition here for speedup
	for (auto e_idx : adj_edges) {

	  Edge e = this->edges.at(e_idx);
	  std::size_t adj_v_idx = e.get_v1_idx() == current_vertex ? e.get_v2_idx() : e.get_v1_idx();
	  
	  //std::size_t adj_v_idx = el.first;
	  //std::size_t e_idx = el.second;


	  
	  //Edge e = this->edges.at(e_idx);

	  // std::cout << "to: "<< adj_v_idx  << " e_idx: " << e.get_color() << " e_idx: " << e_idx << std::endl;
	  
	  //std::size_t adj_v_idx = e.get_v1_idx() == current_vertex ? e.get_v2_idx() : e.get_v1_idx();
	  //std::size_t a = adj.v_idx;

	  //if (a < current_vertex) { continue; }

	  if (seen.find(adj_v_idx) == seen.end()) {
		t.add_tree_edge(current_vertex, adj_v_idx, e.get_color());
		visited.push(adj_v_idx);
		f = true;
		break;
	  }
	  else if (
		!t.is_root(current_vertex) &&
		t.get_parent(current_vertex) != adj_v_idx &&
		!t.has_child(current_vertex, adj_v_idx) &&
		!t.has_ibe(current_vertex, adj_v_idx) &&
		!t.has_obe(current_vertex, adj_v_idx)
	  ) {
		// TODO: why the has child and not parent test?
		//std::cout << "adding back edge: " << current_vertex << " -> " << a << std::endl;
		t.add_be(current_vertex, adj_v_idx, false, e.get_color());
	  }
	}

	if (!f) { visited.pop(); }
  }

  return t;
  }

}; // namespace biedged
