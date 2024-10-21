#include "./graph.hpp"
#include <cstddef>
#include <stack>
#include <unordered_set>
#include <format>
#include <string>
#include <sys/types.h>

namespace povu::graph {
using namespace povu::graph_types;
namespace pgt = povu::graph_types;


/*
  Edge
  ------
 */
Edge::Edge(std::size_t v1_idx, pgt::v_end v1_end , std::size_t v2_idx, pgt::v_end v2_end) : v1_id{v1_idx}, v1_end{v1_end}, v2_id{v2_idx}, v2_end{v2_end} {}
std::size_t Edge::get_v1_idx() const { return this->v1_id; }
pgt::v_end Edge::get_v1_end() const { return this->v1_end; }
std::size_t Edge::get_v2_idx() const { return this->v2_id; }
pgt::v_end Edge::get_v2_end() const { return this->v2_end; }
pgt::side_n_id_t Edge::get_other_vertex(std::size_t v_id) const {
  if (v1_id == v_id) {
    return pgt::side_n_id_t{v2_end, v2_id};
  }
  else {
    return pgt::side_n_id_t{v1_end, v1_id};
  }
}

/*
  Vertex
  ------
 */
Vertex::Vertex(std::size_t v_id) : v_id{v_id} {}
std::size_t Vertex::id() const { return v_id; }
void Vertex::add_edge_l(std::size_t e_id) { e_l.insert(e_id); }
void Vertex::add_edge_r(std::size_t e_id) { e_r.insert(e_id); }
const std::set<std::size_t>& Vertex::get_edges_l() const { return e_l; }
const std::set<std::size_t>& Vertex::get_edges_r() const { return e_r; }

/*
  Graph
  -----
 */
Graph::Graph() {}
Graph::Graph(std::size_t v_count, std::size_t e_count) {
  this->vertices.reserve(v_count);
  this->edges.reserve(e_count);
}

void Graph::add_vertex(std::size_t v_id) {
  vertices.push_back(Vertex{v_id});
  this->v_id_to_idx_.insert(v_id, vertices.size() - 1);
}

std::size_t Graph::v_idx_to_id(std::size_t v_idx) const {
  return this->v_id_to_idx_.get_key(v_idx);
}

std::size_t Graph::v_id_to_idx(std::size_t v_id) const {
  return this->v_id_to_idx_.get_value(v_id);
}

std::size_t Graph::size() const { return this->vertices.size(); }
std::size_t Graph::edge_count() const { return this->edges.size(); }
const std::set<pgt::side_n_id_t> &Graph::tips() const { return this->tips_; }
const Vertex& Graph::get_vertex_by_idx(std::size_t v_idx) const { return vertices[v_idx]; }
const Vertex &Graph::get_vertex_by_id(std::size_t v_id) const {
  return vertices[this->v_id_to_idx_.get_value(v_id)];
}

const Edge& Graph::get_edge(std::size_t e_id) const { return edges[e_id]; }

void Graph::add_tip(std::size_t v_id, pgt::v_end end) {
  std::size_t v_idx = this->v_id_to_idx_.get_value(v_id);
  this->tips_.insert(pgt::side_n_id_t{end, v_idx} );
}

void Graph::add_edge(std::size_t v1_id, pgt::v_end v1_end, std::size_t v2_id, pgt::v_end v2_end) {
  std::size_t v1_idx = this->v_id_to_idx_.get_value(v1_id);
  std::size_t v2_idx = this->v_id_to_idx_.get_value(v2_id);
  edges.push_back(Edge{v1_idx, v1_end, v2_idx, v2_end});
  std::size_t e_idx = edges.size() - 1;

  if (v1_end == pgt::v_end::l) {
    this->vertices[v1_idx].add_edge_l(e_idx);
  }
  else {
    this->vertices[v1_idx].add_edge_r(e_idx);
  }

  if (v2_end == pgt::v_end::l) {
    this->vertices[v2_idx].add_edge_l(e_idx);
  }
  else {
    this->vertices[v2_idx].add_edge_r(e_idx);
  }
}

void Graph::summary() const {
  std::cout << "Bidirected Graph: " << std::endl;
  std::cout << "\t" << "vertex count: " << this->size() << std::endl;
  std::cout << "\t" << "edge count: " << this->edge_count() << std::endl;
}

std::vector<povu::graph::Graph> componetize(const povu::graph::Graph& g, const core::config& app_config) {
  std::string fn_name = std::format("[povu::graph_ops::{}]", __func__);
  if (app_config.verbosity() > 4) { std::cerr << fn_name << std::endl; }

  std::unordered_set<std::size_t> visited, explored;
  std::stack<std::size_t> s;

  s.push(0);
  visited.insert(0);

  std::vector<povu::graph::Graph> components;

  std::set<std::size_t> curr_edges;
  std::set<std::size_t> curr_vertices;
  std::set<pgt::side_n_id_t> curr_tips;

  auto update = [&](std::size_t e_idx, std::size_t curr_v, bool& is_vtx_explored) {
    std::size_t adj_v = g.get_edge(e_idx).get_other_vertex(curr_v).v_idx;
    if (adj_v != curr_v && visited.find(adj_v) == visited.end()) {
      s.push(adj_v);
      visited.insert(adj_v);
      is_vtx_explored = false;
    }
  };

  while (!s.empty()) {
    std::size_t v_idx = s.top();
    const povu::graph::Vertex& v = g.get_vertex_by_idx(v_idx);
    bool is_vtx_explored { true };

    std::set<std::size_t> const& e_l = v.get_edges_l();
    std::set<std::size_t> const& e_r = v.get_edges_r();

    for (std::size_t e_idx: e_l) { update(e_idx, v_idx, is_vtx_explored); }
    for (std::size_t e_idx: e_r) { update(e_idx, v_idx, is_vtx_explored); }


    if (is_vtx_explored) {
      //
      curr_vertices.insert(v.id());
      curr_edges.insert(e_l.begin(), e_l.end());
      curr_edges.insert(e_r.begin(), e_r.end());

      if (g.tips().contains({pgt::v_end::l, v_idx})) {
        curr_tips.insert({pgt::v_end::l, v.id()});
      }
      else if (g.tips().contains({pgt::v_end::r, v_idx})) {
        curr_tips.insert({pgt::v_end::r, v.id()});
      }

      //
      explored.insert(v_idx);
      s.pop();
    }

    if (s.empty()) {
      povu::graph::Graph curr_vg(curr_vertices.size(), curr_edges.size());

      for (std::size_t v_id: curr_vertices) {
        curr_vg.add_vertex(v_id);
      }

      for (std::size_t e_id: curr_edges) {
        const povu::graph::Edge &e = g.get_edge(e_id);
        std::size_t v1_idx = e.get_v1_idx();
        pgt::v_end v1_end = e.get_v1_end();
        std::size_t v1_id = g.v_idx_to_id(v1_idx);

        std::size_t v2_idx = e.get_v2_idx();
        pgt::v_end v2_end = e.get_v2_end();
        std::size_t v2_id = g.v_idx_to_id(v2_idx);

        curr_vg.add_edge(v1_id, v1_end, v2_id, v2_end);
      }

      for (auto [s, id] : curr_tips) {
        curr_vg.add_tip(id, s);
      }

     components.push_back(curr_vg);
     curr_vertices.clear();
     curr_edges.clear();
     curr_tips.clear();
     //curr_vg = povu::graph::Graph();

     for (std::size_t v_idx{}; v_idx < g.size(); ++v_idx) {
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


} // namespace povu::graph
