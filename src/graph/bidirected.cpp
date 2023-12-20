#include <queue>
#include <format>
#include <string>
#include <vector>
#include <utility>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <unordered_set>
#include <map>

#include "./bidirected.hpp"
#include "../core/utils.hpp"
#include "../core/utils.hpp"

namespace bidirected {
/*
 * VertexEnd
 * ---------
 */
std::ostream& operator<<(std::ostream& os, const VertexEnd& ve) {
  switch (ve) {
  case VertexEnd::l:
  os << "+";
  break;
  case VertexEnd::r:
  os << "-";
  break;
  }

  return os;
}


VertexEnd complement(VertexEnd s) {
  return s == VertexEnd::l ? VertexEnd::r : VertexEnd::l;
};


/*
 * Side and SideID
 * ----
 */
bool operator<(const side_n_id_t& lhs, const side_n_id_t& rhs) {
  if (lhs.v_idx < rhs.v_idx) {
	return true;
  }
  else if (lhs.v_idx == rhs.v_idx) {
	return lhs.v_end < rhs.v_end;
  }
  else {
	return false;
  }
}

// implement complement method for side_n_id_t
side_n_id_t side_n_id_t::complement() const {
  return {bidirected::complement(this->v_end), this->v_idx};
}


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

side_n_id_t Edge::get_other_vertex(std::size_t vertex_index) const {
  if (this->get_v1_idx() == vertex_index) {
	return side_n_id_t {this->get_v2_end(), this->get_v2_idx() };
  }
  else {
	return side_n_id_t {this->get_v1_end(), this->get_v1_idx() };
  }
}

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

  this->paths = std::vector<PathInfo>();
  this->handle = std::string();
  this->is_reversed_ = false;
}

Vertex::Vertex(const std::string& label): label(label) {
  this->edges_l = std::set<std::size_t>();
  this->edges_r = std::set<std::size_t>();

  this->paths = std::vector<PathInfo>();
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

std::string Vertex::get_rc_label() const{
  return utils::reverse_complement(this->label);
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

const std::vector<PathInfo>& Vertex::get_paths() const {
  return this->paths;
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

std::vector<side_n_id_t> VariationGraph::get_adj_vertices(std::size_t vertex_index, VertexEnd vertex_end) const {
  std::vector<side_n_id_t> adj_vertices;

  if (vertex_end == VertexEnd::l) {
	for (const auto& edge_index : this->vertices[vertex_index].get_edges_l()) {
	  const Edge& e = this->get_edge(edge_index);
	  adj_vertices.push_back(e.get_other_vertex(vertex_index));
	}
  }
  else {
	for (const auto& edge_index : this->vertices[vertex_index].get_edges_r()) {
		  const Edge& e = this->get_edge(edge_index);
	  adj_vertices.push_back(e.get_other_vertex(vertex_index));
	}
  }

  return adj_vertices;
}

std::unordered_set<id_t> VariationGraph::find_graph_start_nodes() const {
  std::unordered_set<id_t> graph_start_nodes;
  for (id_t i = 0; i < this->size(); ++i) {
	if (this->vertices.at(i).get_edges_l().empty()) {
	  graph_start_nodes.insert(i);
	}
  }
  return graph_start_nodes;
}

std::unordered_set<id_t> VariationGraph::find_graph_end_nodes() const {
  std::unordered_set<id_t> graph_end_nodes;
  for (id_t i = 0; i < this->size(); ++i) {
	if (this->vertices.at(i).get_edges_r().empty()) {
	  graph_end_nodes.insert(i);
	}
  }
  return graph_end_nodes;
}


/**
 * @brief Simplify paths by removing redundant nodes and instead representing it as a single node and strand through which the path passes
 *
 * @param all_paths a vector of paths
 * @return std::vector<subpaths_t> a vector of simplified paths
 */
std::vector<std::vector<side_n_id_t>>
simplify_paths(std::size_t start_id,
			   std::size_t stop_id,
			   const std::vector<std::vector<side_n_id_t>>& all_paths) {

  std::string fn_name = "[povu::graph::VariationGraph::simplify_paths]";
  if (false) {
	std::cerr << fn_name << "\n";
  }

  std::vector<std::vector<side_n_id_t>> simplified_paths;
  simplified_paths.reserve(all_paths.size());

  auto in_flubble = [start_id, stop_id](id_t i) -> bool {
	return i > start_id && i < stop_id;
  };

  for (std::size_t i{}; i < all_paths.size(); ++i) {
	const std::vector<side_n_id_t>& path = all_paths[i];

	if (path.size() % 2 != 0) {
	  std::cout << "Expected even number of vertices in path, skipping path\n";
	  // print path
	  for (auto [side, id]: path) {
		std::cout << side << " " << id << ", ";
	  }
	  std::cout << std::endl;

	  continue;
	}

	std::vector<side_n_id_t> simplified_path;
	simplified_path.reserve(path.size()/2);

	for (std::size_t j{}; j < path.size(); j+=2) {

	  auto [side1, id1] = path[j];
	  auto [side2, id2] = path[j+1];

	  if (side1 == side2) {
		std::cout << "Expected different sides, skipping path\n";
	  }

	  if (id1 != id2) {
		std::cout << "Expected same node id, skipping path\n";
	  }

	  if (!in_flubble(id2)) { continue; }

	  simplified_path.push_back({side1, id2});
	}

	simplified_paths.push_back(simplified_path);
  }

  return simplified_paths;
}


std::vector<std::vector<side_n_id_t>>
VariationGraph::get_paths(id_t start_id, id_t stop_id, bool compact) const {
  std::string fn_name = "[povu::graph::VariationGraph::get_paths]";
  if (false) {
	std::cerr << fn_name << "\n";
  }
  auto in_flubble = [start_id, stop_id](id_t i) -> bool {
	return i >= start_id && i <= stop_id;
  };

  // TODO: remove?
  auto remove = [](std::vector<side_n_id_t>& v, std::function<bool(side_n_id_t)> condition){
	v.erase(std::remove_if(v.begin(), v.end(), condition), v.end());
  };

  // a double ended queue to control the traversal of the region
  std::deque<side_n_id_t> q;

  // a set to keep track of the vertices we've seen
  std::set<side_n_id_t> seen;

  auto unique_push_back_q = [&](side_n_id_t e) {
	if (seen.count(e) || !in_flubble(e.v_idx)) { return; }

	//std::cout << "\t" << "lambda q push back: " << e.v_end << " " << e.v_idx << "\n";

	q.push_back(e);
	seen.insert(e);
  };

  // return the side to which all other vertices connected are greater than i
  // it should be that each side has all vertices greater than i or less than i
  auto greater_side = [&](id_t i) -> std::pair<VertexEnd, int> {
	// vertices adjacent to the left of the vertex at i
	std::vector<side_n_id_t> v_left = this->get_adj_vertices(i, VertexEnd::l);
	// vertices adjacent to the right of the vertex at i
	std::vector<side_n_id_t> v_right = this->get_adj_vertices(i, VertexEnd::r);

	// check that all v_left are less than i
	// or if vertex at i is inverted all v_left are greater than i
	bool l_fwd = std::all_of(v_left.begin(), v_left.end(), [&](side_n_id_t x) { return x.v_idx < i; });
	bool l_rse = std::all_of(v_left.begin(), v_left.end(), [&](side_n_id_t x) { return x.v_idx > i; });


	if (!v_left.empty() && !(l_fwd ^ l_rse)) {
	  return std::make_pair(VertexEnd::l, 1); // TODO: loop boundary
	}

	// check that all v_right are greater than start_idx
	// or if vertex at start_id is inverted all v_right are less then start_idx
	bool r_fwd = std::all_of(v_right.begin(), v_right.end(), [&](side_n_id_t x) { return x.v_idx > i; });
	bool r_rse =std::all_of(v_right.begin(), v_right.end(), [&](side_n_id_t x) { return x.v_idx < i; });

	if (!v_right.empty() && !(r_fwd ^ r_rse)) {
	  return std::make_pair(VertexEnd::l, -2);
	}

	if (l_fwd && r_fwd) { // start vertex is not inverted
	  return std::make_pair(VertexEnd::r, 0);
	}
	else if (l_rse && r_rse) { // start vertex is inverted
	  return std::make_pair(VertexEnd::l, 0);
	}
	else if (l_fwd && r_rse) {
	  return std::make_pair(VertexEnd::r, 2); // TODO: loop boundary
	}
	else {
	  return std::make_pair(VertexEnd::l, -3);
	}
  };

  std::map<side_n_id_t, std::vector<std::vector<side_n_id_t>>> paths_map;
  auto extend_paths = [&paths_map](side_n_id_t frm, side_n_id_t to) {
	paths_map[to] = paths_map[frm];
	for (std::vector<side_n_id_t>& path_ : paths_map[to]) {
	  path_.push_back(to);
	}
  };

  // ==================================

  // ------------------------------------------------
  // Handle the start of the path/flubble/SESE region
  // ------------------------------------------------

  // the start side facing (the rest of) the flubble
  auto [start_side, status_code] = greater_side(start_id);


   
  if (status_code < 0) {
	std::cerr << fn_name << " ERROR: not valid flubble boundary (start) "
			  << start_id << " -> " << stop_id
			  << " status code "  << status_code << "\n";
  }

  const std::vector<side_n_id_t>& to_vertices = this->get_adj_vertices(start_id, start_side);
  // initialize the visit queue with the first elements in the flubble
  std::copy(to_vertices.begin(), to_vertices.end(), std::back_inserter(q));

  // vertices seen before along the path to the key vertex
  std::map<side_n_id_t, std::set<side_n_id_t>> in_path;
  for (std::size_t i{start_id}; i < stop_id+1; ++i) {
	in_path[{VertexEnd::l, i}] = {};
	in_path[{VertexEnd::r, i}] = {};
  }

  in_path.at({start_side, start_id})
	.insert({{v_end_t::l, start_id},
			 {v_end_t::r, start_id}});

  // the key is the vertex id and the value is a vector of (unique) paths
  // leading up to the key vertex


  // start by adding the start vertex as the only path to itself

  paths_map[{start_side, start_id}] =
	std::vector<std::vector<side_n_id_t>>{
	{
	  {complement(start_side), start_id},
	  {start_side, start_id}
	}
  };


  // a pair of side and vertex that has been visited
  std::set<side_n_id_t> visited;
  id_t v_id;
  std::vector<std::vector<side_n_id_t>> curr_paths;
  bidirected::VertexEnd side;

  while (!q.empty()) {

	side_n_id_t side_id_pair = q.front();
	side = side_id_pair.v_end;
	v_id = side_id_pair.v_idx;
	//std::cout << "q.front: (" << side << ", " << v_id << ")\n";

	if (visited.count(side_id_pair)) {
	  q.pop_front();
	  //std::cout << "skipping (already visited)\n";
	  continue;
	}

	//bidirected::Vertex const& v = bi_vg.get_vertex(v_id);

	bool go_back{false};

	{
	  std::vector<bidirected::side_n_id_t> to_adjacents = this->get_adj_vertices(v_id, complement(side));
	  // populate the queue for the next iteration
	  for (auto [s,i]: to_adjacents) {
		unique_push_back_q({s,i});
	  }
	}

	//const std::set<std::size_t>& frm_adj_edges_idxs = get_edges_for_side(v, side);
	  // side == bidirected::VertexEnd::l ? v.get_edges_l() : v.get_edges_r();
	std::vector<side_n_id_t> frm_adjacents;
	for (auto [s, i] : this->get_adj_vertices(v_id, side) ) {
	  frm_adjacents.push_back({s,i});
	}

	//std::cout << "\tfrm adjacents:\n";
	for (side_n_id_t adj: frm_adjacents) {
	  // std::format("\t\tfrom: {}\n", adj.second);
	  //std::cout << "\t\tadj: " << adj.first << ", " << adj.second << "\n";
	  auto [adj_side, adj_id] = adj;
	  //bidirected::VertexEnd adj_side = adj.v_end;
	  v_end_t alt_adj_side = complement(adj_side);


	  if (!in_flubble(adj_id)) { continue; }

	  //std::cout << "\t\t" << adj_side << ", " << alt_adj_side << ", " << adj.second << "\n";
	  side_n_id_t alt_adj = adj.complement();

	  // check whether adj is in paths_map
	  if (paths_map.find(adj) == paths_map.end() && paths_map.find(alt_adj) == paths_map.end() ) {
		// push from to queue and break
		//std::cout << "\t\tpush front: " << adj.first << ", " << adj.second << "\n";
		q.push_front(alt_adj);
		go_back = true;
		break;
	  }

	  // TODO: move this higher
	  // there is a cycle
	  if (in_path[adj].count(side_id_pair)) {
		continue;
	  }

	  if (paths_map.find(adj) == paths_map.end() ) {
		extend_paths(alt_adj, adj);
	  }

	  //subpaths_t& in_paths = paths_map.at(adj);
	  for (std::vector<side_n_id_t>& p: paths_map.at(adj)) {
		std::vector<side_n_id_t> p_ = p;
		p_.push_back(side_id_pair);
		curr_paths.push_back(p_);

		// save that all these vertices are on at least one the paths to v_id
		// populate in_path

		//std::cout << "\t\tinserting:\n";
		for (side_n_id_t x: p_) {
		  //std::cout << "\t\t\t" << x.first << " " << x.second << "\n";
		  in_path[side_id_pair].insert(x);
		}
	  }
	}

	if (!go_back) {
	  paths_map[side_id_pair] = curr_paths;
	  visited.insert(side_id_pair);
	}

	curr_paths.clear();
  }

  auto [stop_side, stop_status_code] = greater_side(stop_id);

  if (stop_status_code < 0) {
	std::cerr << fn_name << " ERROR: not valid flubble boundary (stop) "
			  << start_id << " -> " << stop_id
			  << " status code "  << stop_status_code << "\n";
  }

  //if (stop_status_code == 2) {
  //std::cout << fn_name << stop_side << stop_status_code<< "\n";
  //}


  //std::cerr << fn_name << " exit: " << stop_side << stop_id << "\n";

  side_n_id_t exit = {stop_side, stop_id};

  // if no path has reached the exit, that is because the opposite side
  if (paths_map[exit].empty()) {
	extend_paths(exit.complement(), exit);
  }

  if (false) {
  //std::cout << fn_name << " paths: " << complement(stop_side) << stop_id << "\n";
	for (auto x: paths_map.at(exit) ) {
	  std::cout << fn_name << " path: ";
	  for (auto [s, i]: x) {
		std::cout <<  s << i << ", ";
	  }
	}
	std::cout << "\n";
  }

  if (compact) {
	return simplify_paths(start_id, stop_id, paths_map[exit]);
  }

  return paths_map[exit];
}

const std::vector<path_t>& VariationGraph::get_paths() const {
  return this->paths;
}

void VariationGraph::dbg_print() {
  std::cerr << "VariationGraph: " << std::endl;
  std::cerr << "\t" << "vertex count:  " << this->size() << std::endl;
  std::cerr << "\t" << "valid vertices: " << std::endl;
  std::cerr << "\t" << "edge count: " << this->edges.size() << std::endl;
  std::cerr << "\t" << "path count: " << this->paths.size() << std::endl;

  std::cerr << "\t" << "haplotype start nodes: ";
  utils::print_with_comma(this->haplotype_start_nodes);
  std::cerr << "\n";

  std::cerr << "\t" << "haplotype end nodes: ";
  utils::print_with_comma(this->haplotype_end_nodes);
  std::cerr << "\n";


  std::cerr << "\t" << "graph start nodes: ";
  utils::print_with_comma(this->find_graph_start_nodes());
  std::cerr << "\n";

  std::cerr << "\t" << "graph end nodes: ";
  utils::print_with_comma(this->find_graph_end_nodes());
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

void VariationGraph::add_haplotype_start_node(std::size_t node_id) {
  this->haplotype_start_nodes.insert(node_id);
}

void VariationGraph::add_haplotype_stop_node(std::size_t node_id) {
  this->haplotype_end_nodes.insert(node_id);
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
