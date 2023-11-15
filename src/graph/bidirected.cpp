#include <cstddef>
#include <cstring>
#include <format>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_set>
#include <utility>

#include "./bidirected.hpp"
#include "../core/utils.hpp"

namespace bidirected {
/*
 * Edge
 * ----
 */
Edge::Edge() {
  this->v1_idx = std::size_t();
  this->v1_end = VertexEnd::l;
  this->v2_idx = std::size_t();
  this->v2_end = VertexEnd::l;
}

Edge::Edge(std::size_t v1, VertexEnd v1_end, std::size_t v2, VertexEnd v2_end)
  : v1_idx(v1), v1_end(v1_end),
	v2_idx(v2), v2_end(v2_end)
{}

std::size_t Edge::get_v1_idx() const {
  return this->v1_idx;
}

VertexEnd Edge::get_v1_end() const {
  return this->v1_end;
}

std::size_t Edge::get_v2_idx() const {
  return this->v2_idx;
}

VertexEnd Edge::get_v2_end() const {
  return this->v2_end;
}

// << operator
std::ostream& operator<<(std::ostream& os, const Edge& edge) {
  os << std::format("{{bidirected::Edge {}{} {}{} }}",
					edge.v1_idx, (edge.v1_end == VertexEnd::l ? "+" : "-"),
					edge.v2_idx, (edge.v2_end == VertexEnd::l ? "+" : "-"));

  return os;
}


/*
 * Vertex
 * ------
 */

Vertex::Vertex() {
  this->label = std::string();
  this->edges_l = std::set<std::size_t>();
  this->edges_r = std::set<std::size_t>();

  //std::unordered_set<PathInfo> myset;
  this->paths = std::vector<PathInfo>();
  this->handle = std::string();
  this->is_reversed_ = false;
}

Vertex::Vertex(const std::string& label): label(label) {
  this->edges_l = std::set<std::size_t>();
  this->edges_r = std::set<std::size_t>();

  this->paths = std::vector<PathInfo>();
  //this->paths = std::unordered_set<std::size_t>();
  this->handle = std::string();
  this->is_reversed_ = false;
}

Vertex::Vertex(const std::string& label, const handlegraph::nid_t& id)
  : label(label),
	handle(std::to_string(id))
{
  this->edges_l = std::set<std::size_t>();
  this->edges_r = std::set<std::size_t>();

  this->paths = std::vector<PathInfo>();
  this->is_reversed_ = false;
}

const std::string& Vertex::get_label() const {
  return this->label;
}

const std::string& Vertex::get_handle() const {
  return this->handle;
}

const std::set<std::size_t>& Vertex::get_edges_l() const {
	return this->edges_l;
}

const std::set<std::size_t>& Vertex::get_edges_r() const {
	return this->edges_r;
}

bool Vertex::is_reversed() const {
  return this->is_reversed_;
}

bool Vertex::toggle_reversed() {
  this->is_reversed_ = !this->is_reversed_;
  return this->is_reversed_;
}

void Vertex::add_edge(std::size_t edge_index, VertexEnd vertex_end) {
  if (vertex_end == VertexEnd::l) {
	this->edges_l.insert(edge_index);
  }
  else {
	this->edges_r.insert(edge_index);
  }
}

void Vertex::add_path(std::size_t path_id, std::size_t step_index) {
  this->paths.push_back(PathInfo(path_id, step_index));
}

/*
 * Variation Graph
 * ---------------
 */

VariationGraph::VariationGraph()
  : vertices(std::vector<Vertex>{}),
	edges((std::vector<Edge>{})),
	paths(std::vector<path_t>{})
{}

VariationGraph::VariationGraph(
  std::size_t vertex_count, std::size_t edge_count, std::size_t path_count)
  : vertices(std::vector<Vertex>{}),
	edges(std::vector<Edge>{}),
	paths(std::vector<path_t>{})
{
  this->vertices.reserve(vertex_count);
  this->edges.reserve(edge_count);
  this->paths.reserve(path_count);
}

// getters
// -------
std::size_t VariationGraph::size() const {
  // TODO: should this use some counter instead?
  return this->vertices.size();
}

const Vertex& VariationGraph::get_vertex(std::size_t index) const {
	return this->vertices[index];
}

Vertex& VariationGraph::get_vertex_mut(std::size_t index) {
  return this->vertices.at(index);
}

const Edge& VariationGraph::get_edge(std::size_t index) const {
	return this->edges[index];
}

std::unordered_set<std::size_t> VariationGraph::get_start_nodes() const {
  return this->start_nodes;
}

std::unordered_set<std::size_t> VariationGraph::get_end_nodes() const {
  return this->end_nodes;
}

void VariationGraph::dbg_print() {
  std::cerr << "VariationGraph: " << std::endl;
  std::cerr << "\t" << "vertices: " << this->size() << std::endl;
  std::cerr << "\t"<< "valid vertices: " << std::endl;
  std::cerr << "\t"<< "edges: " << this->edges.size() << std::endl;
  std::cerr << "\t" << "paths: " << this->paths.size() << std::endl;

  std::cerr << "\t" << "start nodes: ";
  utils::print_with_comma(this->start_nodes);
  std::cerr << "\n";

  std::cerr << "\t" << "end nodes: ";
  utils::print_with_comma(this->end_nodes);
  std::cerr << "\n";
}

// setters
// -------
void VariationGraph::append_vertex() {
  this->vertices.push_back(Vertex());
}

void VariationGraph::add_vertex(const Vertex& vertex) {
  this->vertices.push_back(vertex);
}

void VariationGraph::add_edge(std::size_t v1, VertexEnd v1_end,
							  std::size_t v2, VertexEnd v2_end)
{
  std::size_t edge_idx = this->edges.size();
  this->edges.push_back(Edge(v1, v1_end, v2, v2_end));
  this->get_vertex_mut(v1).add_edge(edge_idx, v1_end);
  this->get_vertex_mut(v2).add_edge(edge_idx, v2_end);
}

void VariationGraph::set_min_id(std::size_t min_id) {
  this->min_id = min_id;
}

void VariationGraph::set_max_id(std::size_t max_id) {
  this->max_id = max_id;
}

// HandleGraph
// -----------

bool VariationGraph::has_node(handlegraph::nid_t node_id) const {
  return node_id < this->size();
}

handlegraph::handle_t
VariationGraph::get_handle(const handlegraph::nid_t& node_id, bool is_reverse) const {
  handlegraph::handle_t h;

  std::snprintf(h.data, sizeof(h.data), "%lld", node_id);

  return h;
}

handlegraph::nid_t VariationGraph::get_id(const handlegraph::handle_t& handle) const {
  return std::stoi(handle.data);
}

bool VariationGraph::get_is_reverse(const handlegraph::handle_t& handle) const {
  return this->get_vertex( std::stoll(handle.data)).is_reversed();
}

handlegraph::handle_t VariationGraph::flip(const handlegraph::handle_t& handle) {
  this->get_vertex_mut( std::stoll(handle.data)).toggle_reversed();
  // TODO: handle edge flipping
  return handle;
}

size_t VariationGraph::get_length(const handlegraph::handle_t& handle) const {
  return this->get_vertex(  std::stoll(handle.data)).get_label().length();
}

std::string VariationGraph::get_sequence(const handlegraph::handle_t& handle) const {
  return this->get_vertex(std::stoll(handle.data)).get_label();
}

std::size_t VariationGraph::get_node_count() const {
	return this->size();
}

void VariationGraph::add_start_node(std::size_t node_id) {
  this->start_nodes.insert(node_id);
}

void VariationGraph::add_stop_node(std::size_t node_id) {
  this->end_nodes.insert(node_id);
}

handlegraph::nid_t VariationGraph::min_node_id() const {
	return this->min_id;
}

handlegraph::nid_t VariationGraph::max_node_id() const {
	return this->max_id;
}

bool VariationGraph::follow_edges_impl(
  const handlegraph::handle_t& handle,
  bool go_left,
  const std::function<bool(const handlegraph::handle_t&)>& iteratee) const
{
  return false;
}


bool VariationGraph::for_each_handle_impl(
  const std::function<bool(const handlegraph::handle_t&)>& iteratee,
  bool parallel) const
{
  return false;
}



// MutableHandleGraph
// ------------------
handlegraph::handle_t
VariationGraph::create_handle(const std::string& sequence) {
  handlegraph::handle_t h;

  std::snprintf(h.data, sizeof(h.data), "%ld", this->size());

  this->add_vertex(Vertex(sequence, this->size()));

  return h;
}
handlegraph::handle_t
VariationGraph::create_handle(const std::string& sequence, const handlegraph::nid_t& id) {
  std::size_t id_ = id;

  if (id < this->size()) {
	throw std::invalid_argument("id is less than size");
  }

  // pad with empty/invalid vertices until we reach the id
  for (std::size_t i = this->size(); i < id_; i++) {
	this->append_vertex();
  }

  return create_handle(sequence);
}


//  MutablePathHandleGraph
// -----------------------

handlegraph::path_handle_t
VariationGraph::create_path_handle(const std::string& name, bool is_circular) {
  handlegraph::path_handle_t h;
  std::size_t path_id = this->paths.size();
  this->paths.push_back(path_t({name, path_id, is_circular}));

  strncpy(h.data, std::to_string(path_id).c_str(), sizeof(h.data));

  return h;
}

handlegraph::path_handle_t
VariationGraph::rename_path(const handlegraph::path_handle_t& path_handle,
							const std::string& new_name) {
  // extract path id from path_handle
  std::size_t path_id = std::stoll(path_handle.data);
  this->paths[path_id].name = new_name;
  return path_handle;
}

}; // namespace bidirected
