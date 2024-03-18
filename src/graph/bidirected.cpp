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

std::ostream& operator<<(std::ostream& os, const side_n_id_t& x) {
  os << x.v_idx << x.v_end;
  return os;
}

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

void Edge::set_v1_idx(std::size_t v1_idx) { this->v1_idx = v1_idx; }
void Edge::set_v2_idx(std::size_t v2_idx) { this->v2_idx = v2_idx; }


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

std::set<std::size_t> const& VariationGraph::tips() const {
  return this->tips_;
}

std::set<std::size_t> VariationGraph::graph_start_nodes(bool strict) const {
  // get the intersection of the start nodes and the tips

  // loop over haplotype_start_nodes_ and  extract the vertex ids
  std::set<std::size_t> hap_starts;
  for (auto [_, id] : this->haplotype_start_nodes_){
    hap_starts.insert(id);
  }

  std::set<std::size_t> s;
  std::set_intersection(hap_starts.begin(), hap_starts.end(),
                        this->tips_.begin(), this->tips_.end(),
                        std::inserter(s, s.begin()));
  return s;
}

std::set<std::size_t> VariationGraph::graph_end_nodes(bool strict) const {
  std::set<std::size_t> hap_ends;
  for (auto [_, id] : this->haplotype_end_nodes_){
    hap_ends.insert(id);
  }

  std::set<std::size_t> hap_starts;
  for (auto [_, id] : this->haplotype_start_nodes_){
    hap_starts.insert(id);
  }

  // get the set union of tips_ and graph_end_nodes_ then remove   haplotype_start_nodes_
  std::set<std::size_t> u;
  std::set_union(hap_ends.begin(), hap_ends.end(),
                this->tips_.begin(), this->tips_.end(),
                std::inserter(u, u.begin()));


  std::set<std::size_t> d;
  std::set_difference(u.begin(), u.end(),
                      hap_starts.begin(), hap_starts.end(),
                      std::inserter(d, d.begin()));
  return d;
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
            [](const side_n_id_t &a, const side_n_id_t &b) {
              return a.v_idx < b.v_idx;
  });
  };

  std::cerr << "VariationGraph: " << std::endl;
  std::cerr << "\t" << "vertex count:  " << this->size() << std::endl;
  std::cerr << "\t" << "valid vertices: " << std::endl;
  std::cerr << "\t" << "edge count: " << this->edges.size() << std::endl;
  std::cerr << "\t" << "path count: " << this->paths.size() << std::endl;

  std::cerr << "\t" << "tips: ";
  for (std::size_t v : this->tips_) {
    std::size_t g_id = std::stoull(this->get_vertex(v).get_name());
    this->get_vertex(v).get_edges_l().empty() ?
      buffer.push_back({v_end_t::l, g_id}) :
      buffer.push_back({v_end_t::r, g_id});
    //buffer.push_back({s, g_id});
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
  utils::print_with_comma(std::cerr, this->graph_start_nodes(), ',');
  std::cerr << "\n";

  std::cerr << "\t" << "graph end nodes: ";
  utils::print_with_comma(std::cerr, this->graph_end_nodes(), ',');
  std::cerr << "\n";
}


void VariationGraph::print_dot() const {
  std::cout << "digraph G {" << std::endl;
  std::cout << "\trankdir=LR;" << std::endl;
  std::cout << "\tnode [shape=record];" << std::endl;

  for (std::size_t i{}; i < this->size(); ++i) {
    std::cout << std::format("\t{} [label=\"{}\"];\n", i, this->get_vertex(i).get_name());
  }

  for (std::size_t i{}; i < this->edges.size(); ++i) {
    const Edge& e = this->get_edge(i);
    std::cout << std::format("\t{} -> {};\n", e.get_v1_idx(), e.get_v2_idx());
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


void scc(const VariationGraph &vg) {
  std::string fn_name = std::format("[povu::bidirected::{}]", __func__);

  // a vector of low link values i.e. the smallest vertex id reachable from the
  // current vertex
  // low_llink[i] is the low link value of vertex i
  // low_link[i] = min(dfs_num[i], min(low_link[j] for j in adjacents[i]))
  std::vector<std::size_t> low_link(vg.size(), core::constants::UNDEFINED_SIZE_T);

  std::unordered_set<id_t> visited, explored;
  std::stack<side_n_id_t> s;
  std::size_t vertex_count{};
  std::size_t pre_visit_counter {};

  std::vector<side_n_id_t> tips;
  for (side_n_id_t x : vg.find_haplotype_start_nodes()) {
    tips.push_back(x);
  }
  std::sort(tips.begin(), tips.end(),
            [&vg](const side_n_id_t &a, const side_n_id_t &b) {
              return std::stoull( vg.get_vertex(a.v_idx).get_name()) <
                std::stoull( vg.get_vertex(b.v_idx).get_name());
  });

  for (std::size_t t: vg.tips()) {
    vg.get_vertex(t).get_edges_l().empty() ?
      tips.push_back({VertexEnd::l, t}) :
      tips.push_back({VertexEnd::r, t});
  }
  for (side_n_id_t x : vg.find_haplotype_end_nodes()) {
    tips.push_back(x);
  }
  for (side_n_id_t t : tips) {
    auto [_, v] = t;
    if (visited.find(v) == visited.end()) {
      scc_from_tip(vg, visited, explored, low_link, vertex_count, pre_visit_counter, t);
    }
  }

  for (std::size_t i{}; i < low_link.size(); ++i) {
    if (low_link[i] == core::constants::UNDEFINED_SIZE_T) {
      std::size_t e_idx_l = *vg.get_vertex(i).get_edges_l().begin();
      std::size_t other = vg.get_edge(e_idx_l).get_other_vertex(i).v_idx;
      side_n_id_t t = low_link[other] == core::constants::UNDEFINED_SIZE_T ?
        side_n_id_t{VertexEnd::l, i} :
        side_n_id_t{VertexEnd::r, i};

      scc_from_tip(vg, visited, explored, low_link, vertex_count, pre_visit_counter, t);
    }
  }

  // scc_map for each low link value, store the vertices with that low link value
  std::map<std::size_t, std::vector<std::size_t>> scc_map;
  for (std::size_t i{}; i < low_link.size(); ++i) {
    /*
    if (low_link[i] == core::constants::UNDEFINED_SIZE_T) {
      std::cerr << std::format("{} ERROR: low link value for vertex {} ({}) is undefined\n",
                               fn_name, i, vg.get_vertex(i).get_name());
      exit(1);
    }
    */

    if (scc_map.find(low_link[i]) == scc_map.end()) {
      scc_map[low_link[i]] = std::vector<std::size_t>{i};
    }
    else {
      scc_map[low_link[i]].push_back(i);
    }
  }
}

// a vector of low link values i.e. the smallest vertex id reachable from the
  // current vertex
  // low_llink[i] is the low link value of vertex i
  // low_link[i] = min(dfs_num[i], min(low_link[j] for j in adjacents[i]))
  //std::vector<std::size_t> low_link(vg.size(), core::constants::UNDEFINED_SIZE_T);
//
  // scc_map for each low link value, store the vertices with that low link value
void scc(const VariationGraph &vg,
         std::vector<std::size_t>& low_link,
         std::map<std::size_t, std::vector<std::size_t>>& scc_map) {
  std::string fn_name = std::format("[povu::bidirected::{}]", __func__);


  std::unordered_set<id_t> visited, explored;
  std::stack<side_n_id_t> s;
  std::size_t vertex_count{};
  std::size_t pre_visit_counter {};

  std::vector<side_n_id_t> tips;
  for (side_n_id_t x : vg.find_haplotype_start_nodes()) {
    tips.push_back(x);
  }
  std::sort(tips.begin(), tips.end(),
            [&vg](const side_n_id_t &a, const side_n_id_t &b) {
              return std::stoull( vg.get_vertex(a.v_idx).get_name()) <
                std::stoull( vg.get_vertex(b.v_idx).get_name());
  });

  for (std::size_t t: vg.tips()) {
    vg.get_vertex(t).get_edges_l().empty() ?
      tips.push_back({VertexEnd::l, t}) :
      tips.push_back({VertexEnd::r, t});
  }
  for (side_n_id_t x : vg.find_haplotype_end_nodes()) {
    tips.push_back(x);
  }
  for (side_n_id_t t : tips) {
    auto [_, v] = t;
    if (visited.find(v) == visited.end()) {
      scc_from_tip(vg, visited, explored, low_link, vertex_count, pre_visit_counter, t);
    }
  }

  for (std::size_t i{}; i < low_link.size(); ++i) {
    if (low_link[i] == core::constants::UNDEFINED_SIZE_T) {
      std::size_t e_idx_l = *vg.get_vertex(i).get_edges_l().begin();
      std::size_t other = vg.get_edge(e_idx_l).get_other_vertex(i).v_idx;
      side_n_id_t t = low_link[other] == core::constants::UNDEFINED_SIZE_T ?
        side_n_id_t{VertexEnd::l, i} :
        side_n_id_t{VertexEnd::r, i};

      scc_from_tip(vg, visited, explored, low_link, vertex_count, pre_visit_counter, t);
    }
  }


  for (std::size_t i{}; i < low_link.size(); ++i) {
    /*
    if (low_link[i] == core::constants::UNDEFINED_SIZE_T) {
      std::cerr << std::format("{} ERROR: low link value for vertex {} ({}) is undefined\n",
                               fn_name, i, vg.get_vertex(i).get_name());
      exit(1);
    }
    */

    if (scc_map.find(low_link[i]) == scc_map.end()) {
      scc_map[low_link[i]] = std::vector<std::size_t>{i};
    }
    else {
      scc_map[low_link[i]].push_back(i);
    }
  }
}


std::size_t dfs(const VariationGraph &vg, const std::vector<std::size_t>& tips) {
  std::string fn_name = std::format("[povu::bidirected::{}]", __func__);

  std::size_t vertex_count {};

  std::unordered_set<id_t> visited, explored;
  std::stack<id_t> s;

  for (auto t: tips) {
    s.push(t);
    visited.insert(t);
  }

  //std::size_t f = *vg.graph_start_nodes().begin();
  //s.push(f);
  //visited.insert(f);

  auto update = [&](id_t e_idx, id_t curr_v, bool &is_vtx_explored) {
    id_t adj_v = vg.get_edge(e_idx).get_other_vertex(curr_v).v_idx;
    if (adj_v != curr_v && visited.find(adj_v) == visited.end()) {
      s.push(adj_v);
      visited.insert(adj_v);
      is_vtx_explored = false;
    }
  };

  while (!s.empty()) {
    id_t v = s.top();
    bool is_vtx_explored{true};

    if (explored.find(v) != explored.end()) {
      std::cerr << fn_name << " ERROR: vertex " << v << " already explored\n";
      s.pop();
      continue;
    }

    std::set<std::size_t> const& e_l = vg.get_vertex(v).get_edges_l();
    std::set<std::size_t> const& e_r = vg.get_vertex(v).get_edges_r();

    for (auto e: e_l) { update(e, v, is_vtx_explored); }
    for (auto e: e_r) { update(e, v, is_vtx_explored); }

    if (is_vtx_explored) {
      explored.insert(v);
      s.pop();
      vertex_count++;
    }
  }

  if (vertex_count != vg.size()) {
    std::cerr << fn_name << " ERROR: vertex_count " << vertex_count << " != vg.size() " << vg.size() << std::endl;
    exit(1);
  }

  //std::cerr << fn_name << " " << vertex_count << " " << vg.size() << std::endl;

  return vertex_count;
}

std::vector<std::size_t> sort(const VariationGraph &vg) {
  std::string fn_name = std::format("[povu::bidirected::{}]", __func__);

  std::set<std::size_t> starts { vg.graph_start_nodes() };
  std::set<std::size_t> ends { vg.graph_end_nodes() };

  // left and right side degrees
  std::vector<std::pair<std::size_t, std::size_t>> deg;
  deg.reserve(vg.size());

  for (id_t i = 0; i < vg.size(); ++i) {
    deg.push_back(std::make_pair(vg.get_vertex(i).get_edges_l().size(), vg.get_vertex(i).get_edges_r().size()));
  }

  std::set<std::size_t> seen; // seen vertices TODO: do we need seen edges?
  std::vector<std::size_t> sort_order(vg.size(), core::constants::UNDEFINED_SIZE_T);

  // a pair of id and sort order
  std::vector<std::pair<id_t, id_t>> sort_v;
  sort_v.reserve(vg.size());

  // the side is the side from which the vertex was/would have been reached
  // the vertex is the vertex id
  std::deque<side_n_id_t> q;
  for (id_t idx : starts) {  q.push_back({ deg[idx].first == 0 ? VertexEnd::l : VertexEnd::r , idx}); }

  id_t counter {};

  while (!q.empty()) {
    auto [side, v] = q.front();
    q.pop_front();

    if (seen.find(v) != seen.end()) { continue; }

    sort_v.push_back(std::make_pair(v, counter));
    sort_order.at(v) = counter;
    ++counter;

    Vertex const& v_ref = vg.get_vertex(v);
    std::set<std::size_t> out_edges = side == VertexEnd::l ? v_ref.get_edges_r() : v_ref.get_edges_l();

    // TODO: what if you reach the same vertex from different sides therefore affecting its degree?
    for (std::size_t e_idx: out_edges) {
      auto [side_, v_] = vg.get_edge(e_idx).get_other_vertex(v);

      if (v_ == v) { continue; } // self loop

      std::size_t d = side_ == VertexEnd::l ? --deg.at(v_).first : --deg.at(v_).second;

      std::set<std::size_t> const& eee = side_ == VertexEnd::l ? vg.get_vertex(v_).get_edges_l() : vg.get_vertex(v_).get_edges_r();

      if (d == 0 && seen.find(v_) == seen.end()) {
        q.push_front({ side_, v_ });
        seen.insert(v);
      }
    }
  }

  if (vg.size() != counter + 1) {
    std::cerr << std::format("{} ERROR: graph size {} not equal to vertices found {}",
                             fn_name, vg.size(), counter + 1)
              << std::endl;
    exit(1);
  }

  return sort_order;
}

std::vector<std::size_t> foo(const VariationGraph &vg,
                              const std::vector<std::size_t>& low_link,
                              const std::map<std::size_t, std::vector<std::size_t>>& scc_map) {
  std::string fn_name = std::format("[povu::bidirected::{}]", __func__);

  std::vector<std::pair<std::size_t, std::size_t>> deg; // left and right side degrees
  deg.reserve(vg.size());

  for (id_t i = 0; i < vg.size(); ++i) {
    std::size_t l = vg.get_vertex(i).get_edges_l().size();
    std::size_t r = vg.get_vertex(i).get_edges_r().size();
    std::pair<std::size_t, std::size_t> p = std::make_pair(l, r);
    deg.push_back(p);
  }

  std::set<std::size_t> seen; // seen vertices TODO: do we need seen edges?
  std::vector<std::size_t> sort_order(vg.size(), core::constants::UNDEFINED_SIZE_T);

  utils::TwoWayMap<std::size_t, std::size_t> sort_order_map; // vertex id to sort order

  // a pair of id and sort order
  std::vector<std::pair<id_t, id_t>> sort_v;
  sort_v.reserve(vg.size());

  std::set<side_n_id_t> starts { vg.find_haplotype_start_nodes() };
  std::set<side_n_id_t> tips;
  for (std::size_t t : vg.tips()) {
    vg.get_vertex(t).get_edges_l().empty() ?
      tips.insert({VertexEnd::l, t}) :
      tips.insert({VertexEnd::r, t});
  }

  // get the intersection of tips and starts and store it in a vector bar
  std::vector<side_n_id_t> bar;
  std::set_intersection(starts.begin(), starts.end(), tips.begin(), tips.end(), std::back_inserter(bar));
  std::sort(bar.begin(), bar.end(),
            [&vg](const side_n_id_t &a, const side_n_id_t &b) {
              return std::stoull( vg.get_vertex(a.v_idx).get_name()) < std::stoull( vg.get_vertex(b.v_idx).get_name());
  });

  // the side is the side from which the vertex was/would have been reached
  // the vertex is the vertex id
  std::deque<side_n_id_t> q;
  for (auto x : bar) {  q.push_back(x); }

  {
    std::cerr << "queue starts ";
    for (auto q_ : q) {
      std::cerr << vg.get_vertex(q_.v_idx).get_name() << q_.v_end << ", ";
    }
    std::cerr << std::endl;
  }

  id_t counter {};

  auto red_deg = [&vg, &deg](const side_n_id_t &v) {
    return v.v_end == VertexEnd::l ?
      (deg.at(v.v_idx).first > 0 ? --deg.at(v.v_idx).first : 0) :
      (deg.at(v.v_idx).second > 0 ? --deg.at(v.v_idx).second : 0);
  };

  while (!q.empty()) {
    auto [side, v] = q.front();
    q.pop_front();

    if (seen.find(v) != seen.end()) { continue; }

    sort_v.push_back(std::make_pair(v, counter));
    sort_order.at(v) = counter;
    sort_order_map.insert(v, counter);
    std::cerr << fn_name << " v: " << v << "(" << vg.get_vertex(v).get_name() << ") sort idx: " << counter << std::endl;
    ++counter;

    Vertex const& v_ref = vg.get_vertex(v);
    std::set<std::size_t> out_edges = side == VertexEnd::l ? v_ref.get_edges_r() : v_ref.get_edges_l();
    std::set<std::size_t> in_edges = side == VertexEnd::l ? v_ref.get_edges_l() : v_ref.get_edges_r();

    // TODO: what if you reach the same vertex from different sides therefore affecting its degree?
    for (std::size_t e_idx: out_edges) {
      side_n_id_t other = vg.get_edge(e_idx).get_other_vertex(v);
      auto [side_, v_] = other;

      if (v_ == v) { continue; } // self loop

      std::size_t d = red_deg(other);

      if (d == 0 && seen.find(v_) == seen.end()) {
        q.push_front({ side_, v_ });
        seen.insert(v);
      }
    }

    if (false) {
    if (q.empty()) {
      for (std::size_t e_idx: out_edges) {
        auto [side_, v_] = vg.get_edge(e_idx).get_other_vertex(v);
        //std::size_t d = side_ == VertexEnd::l ? deg.at(v_).first : deg.at(v_).second;

        if (v_ == v) { continue; } // self loop

        if (seen.find(v_) == seen.end()) {
          q.push_back({ side_, v_ });
          seen.insert(v);
        }
      }
    }

    if (q.empty()) {
      for (std::size_t e_idx: in_edges) {
        side_n_id_t other = vg.get_edge(e_idx).get_other_vertex(v);
        auto [side_, v_] = other;
        // std::size_t d = side_ == VertexEnd::l ? deg.at(v_).first : deg.at(v_).second;

        if (v_ == v) { continue; } // self loop

        red_deg(other);
        if (seen.find(v_) == seen.end()) {
          q.push_back({ side_, v_ });
          seen.insert(v);
        }
      }
    }
  }
  }

  std::size_t last =  sort_order_map.get_key(counter-1);
  std::cerr << fn_name << " last " << last << " (" << vg.get_vertex(last).get_name() << ")" << std::endl;

  if (vg.size() != counter) {
    std::cerr << fn_name
              << " ERROR: graph size " << vg.size()
              << " not equal to vertices found " << counter
              << std::endl;
    exit(1);
  }

  return sort_order;
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

std::map<std::size_t, component> VariationGraph::count_components(const core::config& app_config) {
  std::string fn_name = std::format("[povu::bidirected::{}]", __func__);
  if (app_config.verbosity() > 4) { std::cerr << fn_name << std::endl; }

  std::set<std::size_t> visited, explored;
  std::size_t component_count{};
  std::stack<std::size_t> s;

  s.push(0);
  visited.insert(0);

  std::map<std::size_t, component> component_map;

  std::set<side_n_id_t> const& hap_starts = this->find_haplotype_start_nodes();
  std::set<side_n_id_t> const& hap_ends = this->find_haplotype_end_nodes();

  bidirected::VariationGraph curr_vg;
  std::set<std::size_t> curr_paths;
  std::set<std::size_t> curr_edges;
  std::map<std::size_t, std::size_t> vertex_map; // old to new vertex indexes
  std::vector<std::size_t> tips;

  auto update = [&](std::size_t e_idx, std::size_t curr_v, bool& is_vtx_explored) {
    std::size_t adj_v = this->get_edge(e_idx).get_other_vertex(curr_v).v_idx;
    if (adj_v != curr_v && visited.find(adj_v) == visited.end()) {
      s.push(adj_v);
      visited.insert(adj_v);
      is_vtx_explored = false;
    }
  };

  while (!s.empty()) {
    std::size_t v = s.top();
    bool is_vtx_explored { true };

    std::set<std::size_t> const& e_l = this->get_vertex(v).get_edges_l();
    std::set<std::size_t> const& e_r = this->get_vertex(v).get_edges_r();

    for (std::size_t e_idx: e_l) { update(e_idx, v, is_vtx_explored); }
    for (std::size_t e_idx: e_r) { update(e_idx, v, is_vtx_explored); }

    if (is_vtx_explored) {
      explored.insert(v);
      s.pop();

      for (auto p : this->get_vertex(v).get_paths()) {
        curr_paths.insert(p.path_id);
      }

      std::size_t v_ = curr_vg.add_vertex(this->get_vertex(v));
      Vertex& v_mut = curr_vg.get_vertex_mut(v_);
      v_mut.clear_edges();
      v_mut.set_handle(v_+1);
      if (vertex_map.find(v) == vertex_map.end()) {
        vertex_map[v] = v_;
      } else {
        std::cerr << "ERROR: vertex_map already contains " << v << std::endl;
        exit(1);
      }

      if (component_map.find(component_count) == component_map.end() ) {
        component_map[component_count] = component(component_count);
      }

      {
        if (e_l.empty() || e_r.empty()) { // v is a tip or an isolated vertex
          curr_vg.add_tip(v_);
        }

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

      component_map.at(component_count).vertices.push_back(v);
    }

    if (s.empty()) {
      for (std::size_t e: curr_edges) {
        std::size_t v1_idx = this->get_edge(e).get_v1_idx();
        std::size_t v2_idx = this->get_edge(e).get_v2_idx();

        std::size_t e_idx =
          curr_vg.add_edge(vertex_map[this->get_edge(e).get_v1_idx()], this->get_edge(e).get_v1_end(),
                           vertex_map[this->get_edge(e).get_v2_idx()], this->get_edge(e).get_v2_end());
        auto e_ = e_idx + 1;
      }
     curr_edges.clear();

     // add paths to curr_vg
     for (std::size_t p: curr_paths) { curr_vg.add_path(this->get_path(p)); }
     curr_paths.clear();

     component_map.at(component_count).vg = curr_vg;

     std::cerr << fn_name << " component id: " << component_count << std::endl;
     // bidirected::dfs(curr_vg, tips); // for debugging/validating
     curr_vg.dbg_print();

     std::vector<std::size_t> low_link(curr_vg.size(), core::constants::UNDEFINED_SIZE_T);
     std::map<std::size_t, std::vector<std::size_t>> scc_map;
     scc(curr_vg, low_link, scc_map);

     foo(curr_vg, low_link, scc_map);

     tips.clear();
     ++component_count;
     curr_vg = bidirected::VariationGraph();
      for (std::size_t v_idx{}; v_idx < this->size(); ++v_idx) {
        if (explored.find(v_idx) == explored.end()) {
          s.push(v_idx);
          visited.insert(v_idx);
          break;
        }
      }
    }
  }

  return component_map;
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
        if (e_l.empty() || e_r.empty()) { // v is a tip or an isolated vertex
          curr_vg.add_tip(v_);
        }

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

void VariationGraph::add_tip(std::size_t node_id) {
  this->tips_.insert(node_id);
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
