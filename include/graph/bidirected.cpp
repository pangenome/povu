#include "./bidirected.hpp"
#include "spanning_tree.hpp"
#include <cstddef>
#include <cstdint>
#include <unordered_set>
#include <utility>
#include <vector>



namespace povu::bidirected {
using namespace povu::graph_types;
namespace pgt = povu::graph_types;
namespace pbd = povu::bidirected;


/*
  Edge
  ------
 */
Edge::Edge(std::size_t v1_idx, pgt::v_end v1_end , std::size_t v2_idx, pgt::v_end v2_end)
  : v1_idx_{v1_idx}, v1_end{v1_end}, v2_idx_{v2_idx}, v2_end{v2_end} {}
std::size_t Edge::get_v1_idx() const { return this->v1_idx_; }
pgt::v_end Edge::get_v1_end() const { return this->v1_end; }
std::size_t Edge::get_v2_idx() const { return this->v2_idx_; }
pgt::v_end Edge::get_v2_end() const { return this->v2_end; }
pgt::side_n_id_t Edge::get_other_vertex(std::size_t v_idx) const {
  if (this->v1_idx_ == v_idx) {
    return pgt::side_n_id_t{this->v2_end, this->v2_idx_};
  }
  else {
    return pgt::side_n_id_t{this->v1_end, this->v1_idx_};
  }
}

/*
  Vertex
  ------
 */
Vertex::Vertex(std::size_t v_id, const std::string& label) : v_id{v_id}, label_(label) {}
std::size_t Vertex::id() const { return v_id; }
const std::string& Vertex::get_label() const {return this->label_; }
const std::set<std::size_t>& Vertex::get_edges_l() const { return e_l; }
const std::set<std::size_t>& Vertex::get_edges_r() const { return e_r; }
/* setters */
void Vertex::add_edge_l(std::size_t e_id) { e_l.insert(e_id); }
void Vertex::add_edge_r(std::size_t e_id) { e_r.insert(e_id); }
void Vertex::add_ref(std::size_t path_id, pgt::or_t strand, std::size_t step_index) {
  this->refs_.push_back(PathInfo(path_id, strand, step_index));
}


/*
  Graph
  -----
 */

VariationGraph::VariationGraph(std::size_t v_count, std::size_t e_count) {
  this->vertices.reserve(v_count);
  this->edges.reserve(e_count);
  this->dummy_idx_ = constants::UNDEFINED_IDX;
}

std::size_t VG::add_vertex(std::size_t v_id, const std::string& label) {
  vertices.push_back(Vertex{v_id, label});
  this->v_id_to_idx_.insert(v_id, vertices.size() - 1);
  return vertices.size() - 1;
}

std::size_t VG::v_idx_to_id(std::size_t v_idx) const {
  return this->v_id_to_idx_.get_key(v_idx);
}

std::size_t VG::v_id_to_idx(std::size_t v_id) const {
  return this->v_id_to_idx_.get_value(v_id);
}

std::size_t VG::size() const { return this->vertices.size(); }
std::size_t VG::edge_count() const { return this->edges.size(); }
const std::set<pgt::side_n_id_t> &VG::tips() const {
  return this->tips_;
}
const Vertex& VG::get_vertex_by_idx(std::size_t v_idx) const { return vertices[v_idx]; }
const Vertex &VG::get_vertex_by_id(std::size_t v_id) const {
  return vertices[this->v_id_to_idx_.get_value(v_id)];
}
Vertex& VG::get_vertex_mut_by_id(std::size_t v_id) {
    return vertices[this->v_id_to_idx_.get_value(v_id)];
}

const Vertex &VG::get_dummy_vertex() const {
  if (this->dummy_idx_ == constants::UNDEFINED_IDX) {
    throw std::runtime_error("Graph has no dummy vertex");
  }

  return vertices[this->dummy_idx_];
}

std::size_t VG::get_dummy_idx() const {
  return this->dummy_idx_;
}

bool VG::has_dummy() const {
  return this->dummy_idx_ != constants::UNDEFINED_IDX;
}

const Edge& VG::get_edge(std::size_t e_id) const { return edges[e_id]; }

void VG::add_tip(std::size_t v_id, pgt::v_end end) {
  this->tips_.insert( pgt::side_n_id_t{end, v_id} );
}

void VG::add_edge(std::size_t v1_id, pgt::v_end v1_end, std::size_t v2_id, pgt::v_end v2_end) {
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

void VG::add_ref(const std::string& ref_name) {
  std::size_t ref_id = this->refs_.size();
  this->refs_[ref_id] = ref_name;
}

void VG::untip() {

  // Do not add a second dummy vertex
  if (this->tips_.empty() || this->dummy_idx_ != constants::UNDEFINED_IDX) {
    return;
  }

  this->dummy_idx_ = this->add_vertex(pc::DUMMY_VTX_ID, "");

  for (auto [s, v_id]: this->tips_) {
    this->add_edge(pc::DUMMY_VTX_ID, pgt::v_end::r, v_id, s);
  }
}

void VG::shrink_to_fit() {
  this->vertices.shrink_to_fit();
  this->edges.shrink_to_fit();
}

void VG::summary() const {
  std::cout << "Bidirected Graph: " << std::endl;
  std::cout << "\t" << "vertex count: " << this->size() << std::endl;
  std::cout << "\t" << "edge count: " << this->edge_count() << std::endl;
}

void VG::print_dot(std::ostream& os) const {

  // map v end left and right to dot west and east for rectangular vertices
  auto v_end_to_dot = [](pgt::v_end e) {
    return e == pgt::v_end::l ? "w" : "e";
  };

  os << "graph G {" << std::endl;
  os << "\t" << "graph [rankdir=LR];" << std::endl;
  os << "\t" << "node [shape=rectangle];" << std::endl;
  for (size_t v_idx {}; v_idx < this->size(); ++v_idx) {
    const Vertex& v = this->get_vertex_by_idx(v_idx);
    std::string v_id = v.id() == constants::UNDEFINED_ID ? "d" : std::to_string(v.id());

    os << "\t" << v_idx << std::format(" [label=\"+ {} -\"];", v_id )<< std::endl;
  }

  for (const Edge& e: this->edges) {
    os << "\t"
       << e.get_v1_idx() << ":" << v_end_to_dot(e.get_v1_end())
       << " -- "
       << e.get_v2_idx() << ":" << v_end_to_dot(e.get_v2_end()) << ";"
       << std::endl;
  }

  os << "}" << std::endl;
}


// TODO DRY this up
std::vector<VG *> componetize(const povu::bidirected::VG &g) {
  std::string fn_name = std::format("[povu::graph_ops::{}]", __func__);
  //if (app_config.verbosity() > 4) { std::cerr << fn_name << std::endl; }

  std::set<pt::idx_t> visited;
  std::stack<pt::idx_t> s;

  pt::idx_t start_vtx { 0 };

  // avoids creating multiple edges between the same vertices
  std::unordered_set<pt::idx_t> added_edges;

  s.push(start_vtx);
  visited.insert(start_vtx);

  std::vector<VG *> components;
  VG *curr_vg;

  bool is_vtx_explored { true };

  while (!s.empty()) {
    pt::idx_t v_idx = s.top();

    const Vertex& v = g.get_vertex_by_idx(v_idx);

    is_vtx_explored = true;

    for (auto e_idx : v.get_edges_l()) {
      const Edge& e = g.get_edge(e_idx);
      auto [_, adj_v_idx] = e.get_other_vertex(v_idx);
      if (adj_v_idx != v_idx && visited.find(adj_v_idx) == visited.end()) {
        s.push(adj_v_idx);
        visited.insert(adj_v_idx);
        is_vtx_explored = false;
      }
    }

    for (auto e_idx : v.get_edges_r()) {
      const Edge& e = g.get_edge(e_idx);
      auto [_, adj_v_idx] = e.get_other_vertex(v_idx);
      if (adj_v_idx != v_idx && visited.find(adj_v_idx) == visited.end()) {
        s.push(adj_v_idx);
        visited.insert(adj_v_idx);
        is_vtx_explored = false;
      }
    }

    if (is_vtx_explored) {
      s.pop();
    }

    if (s.empty()) {

      curr_vg = new VG(visited.size(), added_edges.size());

      /* add vertices */
      for (auto v_idx : visited) {
        const Vertex& v = g.get_vertex_by_idx(v_idx);
        curr_vg->add_vertex(v.id(), v.get_label());
      }

      /* add edges */
      for (auto v_idx : visited) {
        const Vertex& v = g.get_vertex_by_idx(v_idx);
        for (auto e_idx : v.get_edges_l()) {
          if (added_edges.contains(e_idx)) {
            continue;
          }

          added_edges.insert(e_idx);

          const Edge& e = g.get_edge(e_idx);
          auto [side, adj_v_idx] = e.get_other_vertex(v_idx);
          curr_vg->add_edge(v.id(), pgt::v_end::l, g.v_idx_to_id(adj_v_idx), side);
        }

        for (auto e_idx : v.get_edges_r()) {
          if (added_edges.contains(e_idx)) {
            continue;
          }

          added_edges.insert(e_idx);
          const Edge& e = g.get_edge(e_idx);
          auto [side, adj_v_idx] = e.get_other_vertex(v_idx);
          curr_vg->add_edge(v.id(), pgt::v_end::r, g.v_idx_to_id(adj_v_idx), side);
        }
      }

      /* add tips */
      for (auto [side, v_id] : g.tips()) {
        if ( visited.contains(g.v_id_to_idx(v_id))) {
          curr_vg->add_tip(v_id, side);
        }
      }

      // clear the set for the next component
      added_edges.clear();
      components.push_back(curr_vg);

      /* find the next unvisited vertex */
      for (std::size_t v_idx{}; v_idx < g.size(); ++v_idx) {
        if (!visited.contains(v_idx)) { // if not visited
          visited.clear();
          s.push(v_idx);
          visited.insert(v_idx);
          break;
        }
      }
    }
  }

  return components;
}

pst::Tree compute_spanning_tree(const VG &g) {

  /* assumes that the graph has a dummy vertex */
  if (g.get_dummy_idx() == constants::UNDEFINED_IDX) {
    std::cerr
      << "Graph has no dummy vertex computation of spanning tree expects dummy"
      << std::endl;
    pst::Tree t{ 0 };
  }

  // we bi-edge all vertices except for the single dummy start vertex
  pt::idx_t tree_size = ((g.size() - 1) * 2) + 1;
  pst::Tree t { tree_size };

  std::stack<side_n_id_t> s; // stores the vertex index not id
  std::set<pgt::side_n_id_t> explored;
  std::set<pgt::side_n_id_t> visited;

  std::map<side_n_id_t, pt::idx_t> v_idx_to_tree_idx;

  side_n_id_t start = { pgt::v_end::r, pc::DUMMY_VTX_ID };

  std::set<std::pair<id_t, id_t>> added_edges;

  auto mark_edge = [&](side_n_id_t a, side_n_id_t b) {
    added_edges.insert({v_idx_to_tree_idx[a], v_idx_to_tree_idx[b]});
    added_edges.insert({v_idx_to_tree_idx[b], v_idx_to_tree_idx[a]});
  };

  auto is_marked = [&](side_n_id_t a, side_n_id_t b) -> bool {
    return added_edges.contains(
        std::make_pair(v_idx_to_tree_idx[a], v_idx_to_tree_idx[b]));
  };

  s.push(start);
  visited.insert(start);
  // workaround for dummy vertex being bidirected
  visited.insert( {pgt::v_end::l, pc::DUMMY_VTX_ID} );
  explored.insert( {pgt::v_end::l, pc::DUMMY_VTX_ID} );

  bool is_vtx_explored { true };

  auto get_alt = [](pgt::v_end s) {
    return s == pgt::v_end::l ? pgt::v_end::r : pgt::v_end::l;
  };

  auto end2typ = [](pst::v_end e) {
    return (e == pst::v_end::l) ? v_type_e::l : v_type_e::r;
  };

  pt::idx_t counter { 0 };

  while (!s.empty()) {
    side_n_id_t curr = s.top();
    auto [syd, v_id] = curr;
    pt::idx_t v_idx = g.v_id_to_idx(v_id);

    std::cerr << "curr: " << curr << std::endl;

    const Vertex &v = g.get_vertex_by_id(v_id);

    is_vtx_explored = true;

    if (syd == pgt::v_end::l) {
      for (auto e_idx : v.get_edges_l()) {
        const Edge &e = g.get_edge(e_idx);
        auto [n_adj_side, n_v_idx] = e.get_other_vertex(v_idx);

        pt::id_t n_v_id = g.v_idx_to_id(n_v_idx);
        side_n_id_t neigbour = {n_adj_side, n_v_id};

        v_end_t n_alt_side = get_alt(n_adj_side);
        side_n_id_t n_alt = {n_alt_side, n_v_id};

        if ( !visited.contains(neigbour) ) {

          s.push(neigbour);
          s.push(n_alt);

          if (!v_idx_to_tree_idx.contains(curr)) {
            t.add_vertex({counter++, v_id, end2typ(syd)});
            v_idx_to_tree_idx[curr] = counter - 1;
          }

          t.add_vertex({counter++, n_v_id, end2typ(n_adj_side)});
          v_idx_to_tree_idx[neigbour] = counter - 1;

          t.add_vertex({counter++, n_v_id, end2typ(n_alt_side)});
          v_idx_to_tree_idx[n_alt] = counter - 1;

          t.add_tree_edge(v_idx_to_tree_idx[curr],
                          v_idx_to_tree_idx[neigbour],
                          color::gray);

          

          t.add_tree_edge(v_idx_to_tree_idx[neigbour],
                          v_idx_to_tree_idx[n_alt],
                          color::black);

          mark_edge(curr, neigbour);
          mark_edge(neigbour, n_alt);

          visited.insert(neigbour);
          visited.insert(n_alt);

          is_vtx_explored = false;
          break;
        }
        else {
          // add a backedge if:
          //  - not a parent child
          //  - does not exist

          if (!is_marked(curr, neigbour)) {
            t.add_be(v_idx_to_tree_idx[curr], v_idx_to_tree_idx[neigbour],
                     pst::EdgeType::back_edge, color::gray);
            mark_edge(curr, neigbour);
          }
          
        }
      }
    }

    if (syd == pgt::v_end::r) {
      for (auto e_idx : v.get_edges_r()) {
        const Edge &e = g.get_edge(e_idx);
        auto [n_adj_side, n_v_idx] = e.get_other_vertex(v_idx);

        pt::id_t n_v_id = g.v_idx_to_id(n_v_idx);
        side_n_id_t neigbour = {n_adj_side, n_v_id};

        v_end_t n_alt_side = get_alt(n_adj_side);
        side_n_id_t n_alt = { n_alt_side, n_v_id };

        if ( !visited.contains(neigbour) ) {

          s.push(neigbour);
          s.push(n_alt);

          if (!v_idx_to_tree_idx.contains(curr)) {
            t.add_vertex({counter++, v_id, end2typ(syd)});
            v_idx_to_tree_idx[curr] = counter - 1;
          }

          t.add_vertex({counter++, n_v_id, end2typ(n_adj_side)});
          v_idx_to_tree_idx[neigbour] = counter - 1;

          t.add_vertex({counter++, n_v_id, end2typ(n_alt_side)});
          v_idx_to_tree_idx[n_alt] = counter - 1;

          t.add_tree_edge(v_idx_to_tree_idx[curr],
                          v_idx_to_tree_idx[neigbour],
                          color::gray);

          t.add_tree_edge(v_idx_to_tree_idx[neigbour],
                          v_idx_to_tree_idx[n_alt],
                          color::black);

          mark_edge(curr, neigbour);
          mark_edge(neigbour, n_alt);

          visited.insert(neigbour);
          visited.insert(n_alt);

          is_vtx_explored = false;
          break;
        }
        else {
          // add a backedge if:
          //  - not a parent child
          //  - does not exist

          if (!is_marked(curr, neigbour)) {
            t.add_be(v_idx_to_tree_idx[curr], v_idx_to_tree_idx[neigbour],
                     pst::EdgeType::back_edge, color::gray);
            mark_edge(curr, neigbour);
          }
        }
      }
    }

    if (is_vtx_explored) {
      explored.insert(curr);
      s.pop();
    }
  }

  return t;
}

} // namespace povu::bidirected
