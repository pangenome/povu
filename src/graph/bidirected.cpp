#include <cstddef>
#include <cstring>
#include <format>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <ostream>
#include <stack>
#include <deque>
#include <string>
#include <sys/types.h>
#include <unordered_set>
#include <utility>
#include <vector>


#include "./bidirected.hpp"
#include "../core/utils.hpp"
#include "../core/utils.hpp"


namespace bidirected {

// Definition of operator<< outside the struct
std::ostream& operator<<(std::ostream& os, const component& comp) {
  os << "Component index: " << comp.idx << std::endl;

  os << "Starts: ";
  utils::print_with_comma(os, comp.starts, ',');
  os << std::endl;

  os << "Ends: ";
  utils::print_with_comma(os, comp.ends, ',');
  os << std::endl;

  os << "Orphan tips: ";
  utils::print_with_comma(os, comp.orphan_tips, ',');
  os << std::endl;

  os << "Non-tip haplotype starts: ";
  utils::print_with_comma(os, comp.non_tip_hap_starts, ',');
  os << std::endl;

  os << "Vertices: ";
  utils::print_with_comma(os, comp.vertices, ',');
  os << std::endl;
  return os;
}

/*
 * Orientation
 * ------------
 */
// >> and << might be better than + and -
std::ostream& operator<<(std::ostream& os, const orientation_t& o) {
  switch (o) {
  case orientation_t::forward:
  os << "+";
  break;
  case orientation_t::reverse:
  os << "-";
  break;
  }

  return os;
}

/*
 * id and orientation
 * -----------------
 */
std::ostream& operator<<(std::ostream& os, const id_n_orientation_t& x) {
  os << x.v_idx << x.orientation;
  return os;
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

std::size_t Edge::get_eq_class() const { return this->eq_class; }

void Edge::set_v1_idx(std::size_t v1_idx) { this->v1_idx = v1_idx; }
void Edge::set_v2_idx(std::size_t v2_idx) { this->v2_idx = v2_idx; }
void Edge::set_eq_class(std::size_t eq_class) { this->eq_class = eq_class; }

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
  : label(label), handle(std::to_string(id)) {
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

const std::string& Vertex::get_name() const {
  return this->name_;
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

std::size_t Vertex::get_eq_class() const { return this->eq_class; }

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

void Vertex::clear_edges() {
  this->edges_l.clear();
  this->edges_r.clear();
}

void Vertex::add_path(std::size_t path_id, std::size_t step_index) {
  this->paths.push_back(PathInfo(path_id, step_index));
}

void Vertex::set_handle(const std::string& handle) {
  this->handle = handle;
}
void Vertex::set_handle(id_t id) {
  this->handle = std::to_string(id);
}

void Vertex::set_name(const std::string& name) {
  this->name_ = name;
}

void Vertex::set_eq_class(std::size_t eq_class) { this->eq_class = eq_class; }

/*
 * Variation Graph
 * ---------------
 */
VariationGraph::VariationGraph()
    : vertices(std::vector<Vertex>{}),
      edges((std::vector<Edge>{})),
      paths(std::map<id_t, path_t>{})
{}

// TODO: to remove path_count?
VariationGraph::VariationGraph(std::size_t vertex_count, std::size_t edge_count,
                               std::size_t path_count)
    : vertices(std::vector<Vertex>{}),
      edges(std::vector<Edge>{}),
      paths(std::map<id_t, path_t>{})
{
  this->vertices.reserve(vertex_count);
  this->edges.reserve(edge_count);
}

// Getters
// -------
std::size_t VariationGraph::size() const {
  // TODO: should this use some counter instead?
  return this->vertices.size();
}

const Vertex& VariationGraph::get_vertex(std::size_t index) const {
    return this->vertices[index];
}

const Vertex& VariationGraph::get_vertex_by_name(std::string n) const {
  for (const auto& v : this->vertices) {
    if (v.get_name() == n) {
      return v;
    }
  }

  throw std::runtime_error("Vertex not found");
}

std::size_t VariationGraph::get_vertex_idx_by_name(std::string n) const {
  for (std::size_t i = 0; i < this->vertices.size(); i++) {
    if (this->vertices[i].get_name() == n) {
      return i;
    }
  }

  throw std::runtime_error("Vertex not found");
}

Vertex& VariationGraph::get_vertex_mut(std::size_t index) {
  return this->vertices.at(index);
}

const std::vector<Edge>& VariationGraph::get_all_edges() const {
  return this->edges;
}

Edge& VariationGraph::get_edge_mut(std::size_t index) {
  return this->edges.at(index);
}

const Edge& VariationGraph::get_edge(std::size_t index) const {
    return this->edges[index];
}

std::set<side_n_id_t> VariationGraph::get_orphan_tips() const {
  std::set<side_n_id_t> orphan_tips = this->tips();

  for (auto v : this->find_haplotype_start_nodes()) {
    if (orphan_tips.find(v) != orphan_tips.end()) {
      orphan_tips.erase(v);
    }
  }


  for (auto v : this->find_haplotype_end_nodes()) {
    //auto [_, id] = i;
    if (orphan_tips.find(v) != orphan_tips.end()) {
      orphan_tips.erase(v);
    }
  }

  return orphan_tips;
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

std::set<side_n_id_t> const& VariationGraph::tips() const { return this->tips_; }

std::set<side_n_id_t> VariationGraph::graph_start_nodes(bool strict) const {
  // get the intersection of the start nodes and the tips
  std::set<side_n_id_t> gs;
  std::set<side_n_id_t> hs = this->find_haplotype_start_nodes();
  const std::set<side_n_id_t>& t = this->tips();

  // get the intersection of the start nodes and the tips and save it to gs
  std::set_intersection(hs.begin(), hs.end(),
                        t.begin(), t.end(),
                        std::inserter(gs, gs.begin()));

  return gs;
}

std::set<side_n_id_t> VariationGraph::graph_end_nodes(bool strict) const {
  std::set<side_n_id_t> ge;
  std::set<side_n_id_t> he = this->find_haplotype_end_nodes();
  const std::set<side_n_id_t>& t = this->tips();

  std::set_intersection(he.begin(), he.end(),
                        t.begin(), t.end(),
                        std::inserter(ge, ge.begin()));

  return ge;
}

std::set<side_n_id_t> const& VariationGraph::find_haplotype_start_nodes() const {
  return this->haplotype_start_nodes_;
}

std::set<side_n_id_t> const& VariationGraph::find_haplotype_end_nodes() const {
  return this->haplotype_end_nodes_;
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
 if (false) { std::cerr << fn_name << "\n"; }

 std::vector<std::vector<side_n_id_t>> simplified_paths;
 simplified_paths.reserve(all_paths.size());

 auto in_flubble = [start_id, stop_id](id_t i) -> bool {
  return i > start_id && i < stop_id;
 };

 for (std::size_t i{}; i < all_paths.size(); ++i) {
  const std::vector<side_n_id_t>& path = all_paths[i];

 if (path.size() % 2 != 0) {
  std::cerr << "Expected even number of vertices in path, skipping path\n";
      // print path
 for (auto [side, id]: path) {
  std::cerr << side << " " << id << ", ";
 }
 std::cerr << std::endl;

 continue;
 }

 std::vector<side_n_id_t> simplified_path;
 simplified_path.reserve(path.size()/2);

 for (std::size_t j{}; j < path.size(); j+=2) {
  auto [side1, id1] = path[j];
 auto [side2, id2] = path[j+1];

 if (side1 == side2) {
  std::cerr << "Expected different sides, skipping path\n";
 }

 if (id1 != id2) {
  std::cerr << "Expected same node id, skipping path\n";
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
      { complement(start_side), start_id},
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

std::vector<path_t> VariationGraph::get_paths() const {
  std::string fn_name = std::format("[povu::bidirected::{}]", __func__);
  std::vector<path_t> paths_;

  for (auto [_, v]: this->paths) {
    paths_.push_back(v);
  }

  return paths_;
}

const path_t& VariationGraph::get_path(std::size_t path_id) const {
  return this->paths.at(path_id);
}
std::size_t VariationGraph::get_path_count() const {
  return this->paths.size();
}



void VariationGraph::dbg_print() {
  std::vector<side_n_id_t> buffer;

  auto sort_buffer = [&buffer] () {
    std::sort(buffer.begin(), buffer.end(),
              [](const side_n_id_t &a, const side_n_id_t &b) { return a < b; }
              );
  };

  std::cerr << "VariationGraph: " << std::endl;
  std::cerr << "\t" << "vertex count: " << this->size() << std::endl;
  std::cerr << "\t" << "valid vertices: " << std::endl;
  std::cerr << "\t" << "edge count: " << this->edges.size() << std::endl;
  std::cerr << "\t" << "path count: " << this->paths.size() << std::endl;

  std::cerr << "\t" << "tips: ";
  for (auto [s, v]: this->tips()) {
    std::size_t g_id = std::stoull(this->get_vertex(v).get_name());
    buffer.push_back({s, g_id});
  }
  sort_buffer();
  utils::print_with_comma(std::cerr, buffer, ',');
  buffer.clear();
  std::cerr << "\n";

  std::cerr << "\t" << "orphan tips: ";
  for (auto [s, v]: this->get_orphan_tips()) {
    std::size_t g_id = std::stoull(this->get_vertex(v).get_name());
    buffer.push_back({s, g_id});
  }
  sort_buffer();
  utils::print_with_comma(std::cerr, buffer, ',');
  buffer.clear();
  std::cerr << "\n";

  std::cerr << "\t" << "haplotype start nodes: ";
  for (auto [s, v]: this->haplotype_start_nodes_) {
    std::size_t g_id = std::stoull(this->get_vertex(v).get_name());
    buffer.push_back({s, g_id});
  }
  sort_buffer();
  utils::print_with_comma(std::cerr, buffer, ',');
  buffer.clear();
  std::cerr << "\n";

  std::cerr << "\t" << "haplotype end nodes: ";
  for (auto [s, v]: this->haplotype_end_nodes_) {
    std::size_t g_id = std::stoull(this->get_vertex(v).get_name());
    buffer.push_back({s, g_id});
  }
  sort_buffer();
  utils::print_with_comma(std::cerr, buffer, ',');
  buffer.clear();
  std::cerr << "\n";

  std::cerr << "\t" << "graph start nodes: ";
  for (auto [s, v]: this->graph_start_nodes()) {
    std::size_t g_id = std::stoull(this->get_vertex(v).get_name());
    buffer.push_back({s, g_id});
  }
  sort_buffer();
  utils::print_with_comma(std::cerr, buffer, ',');
  buffer.clear();
  std::cerr << "\n";


  std::cerr << "\t" << "graph end nodes: ";
    for (auto [s, v]: this->graph_end_nodes()) {
    std::size_t g_id = std::stoull(this->get_vertex(v).get_name());
    buffer.push_back({s, g_id});
  }
  sort_buffer();
  utils::print_with_comma(std::cerr, buffer, ',');
  buffer.clear();
  std::cerr << "\n";
}


void VariationGraph::print_dot() const {
  std::cout << "digraph G {" << std::endl;
  std::cout << "\trankdir=LR;" << std::endl;
  std::cout << "\tnode [shape=record];" << std::endl;

  for (std::size_t i{}; i < this->size(); ++i) {
    std::cout << std::format("\t{} [label=\"{} {}\"];\n",
                             i, this->get_vertex(i).get_name(),
                             (this->get_vertex(i).get_eq_class() == core::constants::UNDEFINED_SIZE_T
                              ? ""
                              : std::to_string(this->get_vertex(i).get_eq_class() )  ) );
  }

  for (std::size_t i{}; i < this->edges.size(); ++i) {
    const Edge& e = this->get_edge(i);
    std::cout << std::format("\t{} -> {} [label=\"{}\"];\n", e.get_v1_idx(), e.get_v2_idx(),
                             (this->get_edge(i).get_eq_class() == core::constants::UNDEFINED_SIZE_T
                              ? ""
                              : std::to_string(this->get_edge(i).get_eq_class()))
                             );
  }

  std::cout << "}" << std::endl;
}

// setters & modifiers
// -------------------
void VariationGraph::append_vertex() {
  this->vertices.push_back(Vertex());
}

std::size_t VariationGraph::add_vertex(const Vertex& vertex) {
  this->vertices.push_back(vertex);
  return this->vertices.size() - 1;
}

void VariationGraph::add_edge(const Edge& edge) {
  std::size_t edge_idx = this->edges.size();
  this->edges.push_back(edge);
  this->get_vertex_mut(edge.get_v1_idx()).add_edge(edge_idx, edge.get_v1_end());
  this->get_vertex_mut(edge.get_v2_idx()).add_edge(edge_idx, edge.get_v2_end());
}

std::size_t VariationGraph::add_edge(std::size_t v1, VertexEnd v1_end, std::size_t v2, VertexEnd v2_end) {
  std::size_t edge_idx = this->edges.size();
  this->edges.push_back(Edge(v1, v1_end, v2, v2_end));
  this->get_vertex_mut(v1).add_edge(edge_idx, v1_end);
  this->get_vertex_mut(v2).add_edge(edge_idx, v2_end);

  return edge_idx;
}

void VariationGraph::add_path(const path_t& path) {
  std::string fn_name = std::format("[povu::bidirected::{}]", __func__);

  if (this->paths.find(path.id) != this->paths.end()) {
    throw std::invalid_argument(std::format("{} path id {} already exists in the graph.", fn_name, path.id));
  }

  this->paths[path.id] = path_t{path.name, path.id, path.is_circular};
}

void VariationGraph::set_raw_paths(std::vector<std::vector<id_n_orientation_t>> &raw_paths) {
  this->raw_paths = raw_paths;
}

void VariationGraph::set_min_id(std::size_t min_id) {
  this->min_id = min_id;
}

void VariationGraph::set_max_id(std::size_t max_id) {
  this->max_id = max_id;
}



/**
  * Tarjan's algorithm for finding strongly connected components
  *
  * @param vg
  * @return
 */
void scc_from_tip(const VariationGraph &vg,
                  std::unordered_set<id_t> &visited,
                  std::unordered_set<id_t> &explored,
                  std::vector<std::size_t> &low_link,
                  std::size_t &vertex_count,
                  std::size_t &pre_visit_counter,
                  side_n_id_t tip) {
  std::string fn_name = std::format("[povu::bidirected::{}]", __func__);

  std::stack<side_n_id_t> s;
  s.push(tip);
  visited.insert(tip.v_idx);

  while (!s.empty()) {
    auto [side, v] = s.top();
    low_link[v] = pre_visit_counter++;

    bool is_vtx_explored{true};

    if (explored.find(v) != explored.end()) {
      // std::cerr << fn_name << " ERROR: vertex " << v << " already explored\n";
      s.pop();
      continue;
    }

    const std::set<std::size_t> &out_edges =
      side == VertexEnd::l ? vg.get_vertex(v).get_edges_r() : vg.get_vertex(v).get_edges_l();

    for (std::size_t e_idx : out_edges) {
      id_t adj_v = vg.get_edge(e_idx).get_other_vertex(v).v_idx;
      if (adj_v != v && visited.find(adj_v) == visited.end()) {
        s.push({ vg.get_edge(e_idx).get_other_vertex(v).v_end , adj_v});
        visited.insert(adj_v);
        is_vtx_explored = false;
      }
    }

    if (is_vtx_explored) {
      explored.insert(v);
      s.pop();
      vertex_count++;

      for (std::size_t e_idx : out_edges) {
        id_t adj_v = vg.get_edge(e_idx).get_other_vertex(v).v_idx;
        low_link[v] = std::min(low_link[v], low_link[adj_v]);
      }
    }
  }
}

bool VariationGraph::validate_haplotype_paths() {
  for (std::size_t path_idx{}; path_idx < this->raw_paths.size(); path_idx++) {
    std::vector<bidirected::id_n_orientation_t> const &raw_path = this->raw_paths[path_idx];
    for (std::size_t i{}; i < raw_path.size() - 1; i++) {
      bool x{false};
      auto [v1, o1] = raw_path[i];
      auto [v2, o2] = raw_path[i+1];

      std::set<std::size_t> const &s1_edges =
        o1 == bidirected::orientation_t::forward ?
        this->get_vertex(v1).get_edges_r() : this->get_vertex(v1).get_edges_l();

      std::set<std::size_t> const &s2_edges =
        o2 == bidirected::orientation_t::forward ?
        this->get_vertex(v2).get_edges_l() : this->get_vertex(v2).get_edges_r();

      std::set<std::size_t> intersection; // shared edges between v1 and v2
      std::set_intersection(s1_edges.begin(), s1_edges.end(),
                            s2_edges.begin(), s2_edges.end(),
                            std::inserter(intersection, intersection.begin()));

      for (std::size_t e_idx : intersection) {
        Edge const& e = this->get_edge(e_idx);
        bool v_valid =
          (e.get_v1_idx() == v1 && e.get_v2_idx() == v2) ||
          (e.get_v1_idx() == v2 && e.get_v2_idx() == v1);

        VertexEnd s1 = o1 == bidirected::orientation_t::forward ? VertexEnd::r : VertexEnd::l;
        VertexEnd s2 = o2 == bidirected::orientation_t::forward ? VertexEnd::l : VertexEnd::r;

        bool s_valid =
          (e.get_v1_end() == s1 && e.get_v2_end() == s2) ||
          (e.get_v1_end() == s2 && e.get_v2_end() == s1);

        x = x || (s_valid && v_valid);
      }

      if (!x) {
        std::cerr << "ERROR: path " << path_idx << " is not valid at " << raw_path[i] << "," << raw_path[i+1] << std::endl;
        //std::cerr << std::format("ERROR: path {} is not valid at {}, {}", path_idx, raw_path[i], raw_path[i+1] ) << std::endl;
        exit(1);
      }
    }
  }

  return true;
}

std::vector<VariationGraph> componetize(const VariationGraph &vg, const core::config& app_config) {
  std::string fn_name = std::format("[povu::bidirected::{}]", __func__);
  if (app_config.verbosity() > 4) { std::cerr << fn_name << std::endl; }

  std::set<std::size_t> visited, explored;
  std::stack<std::size_t> s;

  s.push(0);
  visited.insert(0);

  std::vector<VariationGraph> components;
  bidirected::VariationGraph curr_vg;

  std::set<side_n_id_t> const& hap_starts = vg.find_haplotype_start_nodes();
  std::set<side_n_id_t> const& hap_ends = vg.find_haplotype_end_nodes();

  std::vector<std::size_t> tips;
  std::set<std::size_t> curr_paths;
  std::set<std::size_t> curr_edges;
  std::map<std::size_t, std::size_t> vertex_map; // old to new vertex indexes

  auto update = [&](std::size_t e_idx, std::size_t curr_v, bool& is_vtx_explored) {
    std::size_t adj_v = vg.get_edge(e_idx).get_other_vertex(curr_v).v_idx;
    if (adj_v != curr_v && visited.find(adj_v) == visited.end()) {
      s.push(adj_v);
      visited.insert(adj_v);
      is_vtx_explored = false;
    }
  };

  while (!s.empty()) {
    std::size_t v = s.top();
    bool is_vtx_explored { true };

    std::set<std::size_t> const& e_l = vg.get_vertex(v).get_edges_l();
    std::set<std::size_t> const& e_r = vg.get_vertex(v).get_edges_r();

    for (std::size_t e_idx: e_l) { update(e_idx, v, is_vtx_explored); }
    for (std::size_t e_idx: e_r) { update(e_idx, v, is_vtx_explored); }

    if (is_vtx_explored) {
      explored.insert(v);
      s.pop();

      for (auto p : vg.get_vertex(v).get_paths()) { curr_paths.insert(p.path_id); }

      std::size_t v_ = curr_vg.add_vertex(vg.get_vertex(v));
      Vertex& v_mut = curr_vg.get_vertex_mut(v_);
      v_mut.clear_edges();
      v_mut.set_handle(v_+1);

      if (vertex_map.find(v) == vertex_map.end()) {
        vertex_map[v] = v_;
      }
      else {
        std::cerr << "ERROR: vertex_map already contains " << v << std::endl;
        exit(1);
      }

      {
        // add tips
        if (e_l.empty() && e_r.empty()) {
          // don't warn because warning was already given when reading the graph
          curr_vg.add_tip(v_, VertexEnd::l);
        }
        else if (e_r.empty()) {
          curr_vg.add_tip(v_, VertexEnd::r);
        }
        else if (e_l.empty()) {
          curr_vg.add_tip(v_, VertexEnd::l);
        }

        //
        if (hap_starts.count({v_end_t::l, v})) {
          curr_vg.add_haplotype_start_node({v_end_t::l, v_});
        }

        if (hap_starts.count({v_end_t::r, v})) {
          curr_vg.add_haplotype_start_node({v_end_t::r, v_});
        }

        if (hap_ends.count({v_end_t::l, v})) {
          curr_vg.add_haplotype_stop_node({v_end_t::l, v_});
        }

        if (hap_ends.count({v_end_t::r, v})) {
          curr_vg.add_haplotype_stop_node({v_end_t::r, v_});
        }

        for (auto e: e_l) { curr_edges.insert(e); }
        for (auto e: e_r) { curr_edges.insert(e); }
      }
    }

    if (s.empty()) {
      for (std::size_t e: curr_edges) {
          curr_vg.add_edge(vertex_map[vg.get_edge(e).get_v1_idx()], vg.get_edge(e).get_v1_end(),
                           vertex_map[vg.get_edge(e).get_v2_idx()], vg.get_edge(e).get_v2_end());
      }

     // add paths to curr_vg
     for (std::size_t p: curr_paths) { curr_vg.add_path(vg.get_path(p)); }
     components.push_back(curr_vg);

     curr_edges.clear();
     curr_paths.clear();
     tips.clear();
     curr_vg = bidirected::VariationGraph();

     for (std::size_t v_idx{}; v_idx < vg.size(); ++v_idx) {
       if (explored.find(v_idx) == explored.end()) {
         s.push(v_idx);
         visited.insert(v_idx);
         break;
       }
     }
    }
  }

  return components;
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

/*
bool VariationGraph::add_graph_start_node(std::size_t node_id) {
  bool l = this->get_vertex(node_id).get_edges_l().empty();
  bool r = this->get_vertex(node_id).get_edges_r().empty();

  if ((l && r) || (!l && !r) ) { return false; }

  VertexEnd v_end = l ? bidirected::VertexEnd::l : bidirected::VertexEnd::r;
  this->graph_start_nodes_.insert(node_id);
  return true;
}

bool VariationGraph::add_graph_end_node(std::size_t node_id) {
  bool l = this->get_vertex(node_id).get_edges_l().empty();
  bool r = this->get_vertex(node_id).get_edges_r().empty();

  if ((l && r) || (!l && !r) ) { return false; }

  VertexEnd v_end = l ? bidirected::VertexEnd::l : bidirected::VertexEnd::r;
  this->graph_end_nodes_.insert(node_id);
  return true;
}
*/
void VariationGraph::add_tip(std::size_t node_id, bidirected::VertexEnd v_end) {
  this->tips_.insert({v_end, node_id});
}

void VariationGraph::add_haplotype_start_node(side_n_id_t i) {
  this->haplotype_start_nodes_.insert(i);
}

void VariationGraph::add_haplotype_stop_node(side_n_id_t i) {
  this->haplotype_end_nodes_.insert(i);
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

handlegraph::path_handle_t VariationGraph::create_path_handle(const std::string& name, bool is_circular) {
  handlegraph::path_handle_t h;
  std::size_t path_id = this->paths.size();
  this->add_path(path_t({name, path_id, is_circular}));

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
