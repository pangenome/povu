#include "biedged.hpp"
#include "bidirected.hpp"
#include <cstddef>
#include <format>
#include <iostream>
#include <string>
#include <vector>

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
	return handle;
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

		bidirected::Edge e_ = g.get_edge(e_idx);

		std::size_t new_v1{};
		VertexType v1_type{};


		if (e_.get_v1_idx() == v_idx){
		  new_v1 = e_.get_v2_end() == bidirected::VertexEnd::l ?
			new_l(e_.get_v2_idx()) :
			new_r(e_.get_v2_idx());

		  v1_type =
			e_.get_v2_end() == bidirected::VertexEnd::l ?
			VertexType::l : VertexType::r;
		} else {
		  new_v1 =
			e_.get_v1_end() == bidirected::VertexEnd::l ?
			new_l(e_.get_v1_idx()) :
			new_r(e_.get_v1_idx());

		  v1_type =
			e_.get_v1_end() == bidirected::VertexEnd::l ?
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

		bidirected::Edge e_ = g.get_edge(e_idx);

		std::size_t new_v2{};
		VertexType v2_type;

		if (e_.get_v2_idx() == v_idx){
		  new_v2 =
			e_.get_v1_end() == bidirected::VertexEnd::l ?
			new_l(e_.get_v1_idx()) :
			new_r(e_.get_v1_idx());

		  v2_type =
		    e_.get_v1_end() == bidirected::VertexEnd::l ?
			VertexType::l : VertexType::r;
		}
		else {
		  new_v2 =
			e_.get_v2_end() == bidirected::VertexEnd::l ?
			new_l(e_.get_v2_idx()) :
			new_r(e_.get_v2_idx());

		  v2_type =
			e_.get_v2_end() == bidirected::VertexEnd::l ? VertexType::l : VertexType::r;
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

}; // namespace biedged
