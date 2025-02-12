#include <cstddef>
#include <cstdint>
#include <iostream>
#include <stack>
#include <string>
#include <sys/types.h>
#include <unordered_set>
#include <utility>
#include <vector>

#include "./bidirected.hpp"


namespace povu::bidirected {

/*
  Edge
  ------
 */
Edge::Edge(pt::idx_t v1_idx, pgt::v_end_e v1_end , pt::idx_t v2_idx, pgt::v_end_e v2_end)
  : v1_idx_{v1_idx}, v1_end{v1_end}, v2_idx_{v2_idx}, v2_end{v2_end} {}
pt::idx_t Edge::get_v1_idx() const { return this->v1_idx_; }
pgt::v_end_e Edge::get_v1_end() const { return this->v1_end; }
pt::idx_t Edge::get_v2_idx() const { return this->v2_idx_; }
pgt::v_end_e Edge::get_v2_end() const { return this->v2_end; }
pgt::side_n_idx_t Edge::get_other_vtx(pt::idx_t v_idx) const {
  return (get_v1_idx() == v_idx) ? pgt::side_n_id_t{v2_end, v2_idx_}
                                 : pgt::side_n_id_t{v1_end, v1_idx_};
}

pgt::side_n_idx_t Edge::get_other_vtx(pt::idx_t v_idx, pgt::v_end_e ve) const {
  const pt::idx_t v1 = get_v1_idx();
  const pt::idx_t v2 = get_v2_idx();

  if (v1 == v2) { // Handle self-loop case
    return {pgt::complement(ve), v1};
  }

  // Return the opposite vertex
  return (v1 == v_idx) ? pgt::side_n_id_t{v2_end, v2_idx_} : pgt::side_n_id_t{v1_end, v1_idx_};
}


/*
  Vertex
  ------
 */
Vertex::Vertex(pt::id_t v_id, const std::string& label) : v_id{v_id}, label_(label) {}
pt::id_t Vertex::id() const { return v_id; }
/* getters */
const std::string &Vertex::get_label() const { return this->label_; }
const std::set<pt::idx_t>& Vertex::get_edges_l() const { return e_l; }
const std::set<pt::idx_t>& Vertex::get_edges_r() const { return e_r; }
/* setters */
void Vertex::add_edge_l(pt::idx_t e_idx) { e_l.insert(e_idx); }
void Vertex::add_edge_r(pt::idx_t e_idx) { e_r.insert(e_idx); }
void Vertex::add_ref(pt::idx_t path_id, pgt::or_t strand, pt::idx_t step_index) {
  this->refs_.push_back(PathInfo(path_id, strand, step_index));
}


/*
  Graph
  -----
 */

VariationGraph::VariationGraph(pt::idx_t v_count, pt::idx_t e_count) {
  this->vertices.reserve(v_count);
  this->edges.reserve(e_count);
}

/* getters */

pt::id_t VG::v_idx_to_id(pt::idx_t v_idx) const {
  return this->v_id_to_idx_.get_key(v_idx);
}

pt::idx_t VG::v_id_to_idx(pt::id_t v_id) const {
  return this->v_id_to_idx_.get_value(v_id);
}

pt::idx_t VG::vtx_count() const { return this->vertices.size(); }

pt::idx_t VG::edge_count() const { return this->edges.size(); }

const std::set<pgt::side_n_id_t> &VG::tips() const {
  return this->tips_;
}
const Edge &VG::get_edge(pt::idx_t e_idx) const { return edges[e_idx]; }
const Vertex& VG::get_vertex_by_idx(pt::idx_t v_idx) const { return vertices[v_idx]; }
const Vertex &VG::get_vertex_by_id(pt::id_t v_id) const {
  return vertices[this->v_id_to_idx_.get_value(v_id)];
}
Vertex& VG::get_vertex_mut_by_id(pt::id_t v_id) {
    return vertices[this->v_id_to_idx_.get_value(v_id)];
}

void VG::add_tip(std::size_t v_id, pgt::v_end_e end) {
  this->tips_.insert( pgt::side_n_id_t{end, v_id} );
}

/* setters */

pt::idx_t VG::add_vertex(pt::id_t v_id, const std::string &label) {
  vertices.push_back(Vertex{v_id, label});
  this->v_id_to_idx_.insert(v_id, vertices.size() - 1);
  return vertices.size() - 1;
}

pt::idx_t VG::add_edge(pt::id_t v1_id, pgt::v_end_e v1_end, pt::id_t v2_id, pgt::v_end_e v2_end) {
  pt::idx_t v1_idx = this->v_id_to_idx_.get_value(v1_id);
  pt::idx_t v2_idx = this->v_id_to_idx_.get_value(v2_id);
  edges.push_back(Edge{v1_idx, v1_end, v2_idx, v2_end});
  pt::idx_t e_idx = edges.size() - 1;

  if (v1_end == pgt::v_end_e::l) {
    this->vertices[v1_idx].add_edge_l(e_idx);
  }
  else {
    this->vertices[v1_idx].add_edge_r(e_idx);
  }

  if (v2_end == pgt::v_end_e::l) {
    this->vertices[v2_idx].add_edge_l(e_idx);
  }
  else {
    this->vertices[v2_idx].add_edge_r(e_idx);
  }

  return e_idx;
}

void VG::add_ref(const std::string& ref_name) {
  std::size_t ref_id = this->refs_.size();
  this->refs_[ref_id] = ref_name;
}

void VG::shrink_to_fit() {
  this->vertices.shrink_to_fit();
  this->edges.shrink_to_fit();
}

void VG::summary() const {
  std::cout << "Bidirected Graph: " << std::endl;
  std::cout << "\t" << "vertex count: " << this->vtx_count() << std::endl;
  std::cout << "\t" << "edge count: " << this->edge_count() << std::endl;
  std::cout << "\t" << "Tip count " << this->tips().size() << std::endl;
  std::cout << "\t";
  pu::print_with_comma(std::cout, this->tips(), ',');
  std::cout << std::endl;
}

void VG::print_dot(std::ostream& os) const {

  // map v end left and right to dot west and east for rectangular vertices
  auto v_end_to_dot = [](pgt::v_end_e e) {
    return e == pgt::v_end_e::l ? "w" : "e";
  };

  os << "graph G {" << std::endl;
  os << "\t" << "graph [rankdir=LR];" << std::endl;
  os << "\t" << "node [shape=rectangle];" << std::endl;
  for (size_t v_idx {}; v_idx < this->vtx_count(); ++v_idx) {
    const Vertex& v = this->get_vertex_by_idx(v_idx);
    std::string v_id = v.id() == constants::UNDEFINED_ID ? "d" : std::to_string(v.id());

    os << "\t" << v_idx
       << std::format(
              "[style=filled, fillcolor=lightblue, label=\"+ {} - \\n ({})\"];",
              v_id, v_idx)
       << std::endl;
  }

  for (const Edge& e: this->edges) {
    os << "\t"
       << e.get_v1_idx() << ":" << v_end_to_dot(e.get_v1_end())
       << " -- "
       << e.get_v2_idx() << ":" << v_end_to_dot(e.get_v2_end())
       << "[color=gray];"
       << std::endl;
  }

  os << "}" << std::endl;
}

std::vector<VG *> componetize(const povu::bidirected::VG &g) {
  std::string fn_name = std::format("[povu::graph_ops::{}]", __func__);

  std::unordered_set<pt::idx_t> visited;
  visited.reserve(g.vtx_count());

  // avoids creating multiple edges between the same vertices
  std::unordered_set<pt::idx_t> added_edges;
  added_edges.reserve(g.edge_count());

  std::stack<pt::idx_t> s;
  pt::idx_t start_vtx{0};
  s.push(start_vtx);
  visited.insert(start_vtx);

  std::vector<VG *> components;
  VG *curr_vg { nullptr };

  std::set<pt::idx_t> comp_vtxs; // current component vertices
  comp_vtxs.insert(start_vtx);

  /* ---------- Helper Functions ---------- */

  auto process_edge = [&](pt::idx_t v_idx, pt::idx_t e_idx) -> void {
    const Edge &e = g.get_edge(e_idx);
    auto [_, adj_v_idx] = e.get_other_vtx(v_idx);

    if (visited.contains(adj_v_idx)) return; // also handles self loops

    s.push(adj_v_idx);
    visited.insert(adj_v_idx);
    comp_vtxs.insert(adj_v_idx);
    return;
  };

  auto add_edges = [&](const Vertex &v, pgt::v_end_e ve, pt::idx_t v_idx, pt::idx_t e_idx) -> void {
    if (added_edges.contains(e_idx)) return; // don't duplicate edges

    added_edges.insert(e_idx);
    const Edge &e = g.get_edge(e_idx);

    auto [adj_s, adj_v_idx] = e.get_other_vtx(v_idx, ve); // handles self loops
    curr_vg->add_edge(v.id(), ve, g.v_idx_to_id(adj_v_idx), adj_s);
  };

  /* ---------- Main Component Search Loop ---------- */

  while (!s.empty()) {
    pt::idx_t v_idx = s.top();
    s.pop();
    const Vertex& v = g.get_vertex_by_idx(v_idx);

    for (auto e_idx : v.get_edges_l()) {
      process_edge(v_idx, e_idx);
    }

    for (auto e_idx : v.get_edges_r()) {
      process_edge(v_idx, e_idx);
    }

    if (s.empty()) {
      curr_vg = new VG(comp_vtxs.size(), added_edges.size());

      /* add vertices */
      for (auto v_idx : comp_vtxs) {
        const Vertex& v = g.get_vertex_by_idx(v_idx);
        curr_vg->add_vertex(v.id(), v.get_label());
      }

      /* add edges */
      for (auto v_idx : comp_vtxs) {
        const Vertex& v = g.get_vertex_by_idx(v_idx);
        for (auto e_idx : v.get_edges_l()) {
          add_edges(v, pgt::v_end_e::l, v_idx, e_idx);
        }

        for (auto e_idx : v.get_edges_r()) {
          add_edges(v, pgt::v_end_e::r, v_idx, e_idx);
        }
      }

      /* add tips */
      for (auto [side, v_id] : g.tips()) {
        if ( comp_vtxs.contains(g.v_id_to_idx(v_id))) {
          curr_vg->add_tip(v_id, side);
        }
      }

      // clear the set for the next component
      components.push_back(curr_vg);
      curr_vg = nullptr;
      added_edges.clear();

      /* find the next unvisited vertex */
      for (std::size_t v_idx{}; v_idx < g.vtx_count(); ++v_idx) {
        if (!visited.contains(v_idx)) { // if not visited
          comp_vtxs.clear();
          s.push(v_idx);
          visited.insert(v_idx);
          comp_vtxs.insert(v_idx);
          break;
        }
      }
    }
  }

  return components;
}

pst::Tree compute_spanning_tree(const VG &g) {

  const bool has_tips{!g.tips().empty()};
  const pt::idx_t root_idx{0};

  std::stack<pt::idx_t> s;
  std::vector<u_int8_t> visited(g.vtx_count(), 0);

  std::set<std::pair<pt::idx_t, pt::idx_t>> added_edges;
  std::unordered_set<pt::idx_t> self_loops;

  pt::idx_t counter{0}; // dfs num
  bool found_new_neighbour { false }; // neighbours exhausted

  pt::idx_t p_idx {pc::INVALID_IDX}; // parent_idx

  pt::idx_t t_vtx_count = has_tips ? (2 * g.vtx_count()) + 1 : 2 * g.vtx_count();
  pst::Tree t{t_vtx_count};

  // biedged idx to tree idx (or counter)
  std::vector<pt::id_t> be_idx_to_ctr(t_vtx_count, 0);

  auto unordered_pair = [](pt::idx_t a, pt::idx_t b) -> std::pair<pt::idx_t, pt::idx_t> {
    return {std::min(a, b), std::max(a, b)};
  };

  auto connect = [&](pt::idx_t a, pt::idx_t b) -> void {
    added_edges.insert(unordered_pair(a, b));
  };

  auto are_connected = [&](pt::idx_t a, pt::idx_t b) -> bool {
    return added_edges.contains(unordered_pair(a, b));
  };

  auto to_be = [&g](pgt::side_n_id_t i) -> pt::idx_t {
    auto [s, v_id] = i;
    pt::idx_t v_idx = g.v_id_to_idx(v_id);

    return (s == pgt::v_end_e::l) ? v_idx * 2 : (v_idx * 2) + 1;
  };

  auto to_bd = [&g](pt::idx_t be_v_idx) -> pgt::side_n_idx_t {
    pgt::v_end_e s = (be_v_idx % 2 == 0) ? pgt::v_end_e::l : pgt::v_end_e::r;
    pt::id_t v_id = g.v_idx_to_id(be_v_idx / 2);
    return {s, v_id};
  };


  auto end2typ = [](pgt::v_end_e e) -> pgt::v_type_e {
    return pgt::v_end_e::l == e ? pgt::v_type_e::l : pgt::v_type_e::r;
  };

  auto add_vertex_to_tree = [&](pgt::v_end_e e, pt::idx_t bd_v_idx) -> void {
    const Vertex &v = g.get_vertex_by_idx(bd_v_idx);

    t.add_vertex({counter++, v.id(), end2typ(e)});
    t.add_vertex({counter++, v.id(), end2typ(pgt::complement(e))});

    be_idx_to_ctr[to_be({e, v.id()})] = counter - 2;
    be_idx_to_ctr[to_be({pgt::complement(e), v.id()})] = counter - 1;

    // add edges
    if (p_idx != pc::INVALID_IDX) {
      t.add_tree_edge(p_idx, counter - 2, pgt::color_e::gray);
      connect(p_idx, counter - 2);
    }
    t.add_tree_edge(counter - 2, counter - 1 , pgt::color_e::black);
    connect(counter-2, counter - 1);
  };

  // returns true if it discovers a new vertex (neighbour), false otherwise
  auto process_edge = [&](pt::idx_t bd_v_idx, pgt::v_end_e ve, pt::idx_t e_idx) -> bool {
    const Edge &e = g.get_edge(e_idx);

    auto [os, ov_idx] = e.get_other_vtx(bd_v_idx, ve); // o for other
    pt::idx_t o_be_idx = to_be( {os, g.v_idx_to_id(ov_idx)});

    if (!visited[ov_idx]) { // has not been visited
      add_vertex_to_tree(os, ov_idx);

      visited[ov_idx] = 1;
      s.push(to_be( {os, g.v_idx_to_id(ov_idx)} ));
      s.push(to_be( {pgt::complement(os), g.v_idx_to_id(ov_idx)} ));

      return true;
    }
    else if (!are_connected(p_idx, be_idx_to_ctr[o_be_idx])) {
      // add a backedge if:
      //  - not a parent child relationship
      //  - a backedge does not already exist
      t.add_be(p_idx, be_idx_to_ctr[o_be_idx], pst::e_type_e::back_edge, pgt::color_e::gray);
      connect(p_idx, be_idx_to_ctr[o_be_idx]);
    }
    else if (__builtin_expect((bd_v_idx == ov_idx && !self_loops.contains(bd_v_idx)), 0)) {
      // add a self loop backedge, a parent-child relationship
      t.add_be(p_idx, be_idx_to_ctr[o_be_idx] , pst::e_type_e::back_edge, pgt::color_e::gray);
      self_loops.insert(bd_v_idx);
    }

    return false;
  };

  if (has_tips) { // add a dummy vertex to the tree
    p_idx = counter;
    t.add_vertex({counter++, pc::DUMMY_VTX_ID, v_type_e::dummy});
  }

  side_n_id_t start = has_tips ? *g.tips().begin() : pgt::side_n_id_t{pgt::v_end_e::l, g.v_idx_to_id(0)};
  auto [s_v_end, s_v_id] = start;
  pt::idx_t s_v_idx = g.v_id_to_idx(s_v_id);
  s.push( to_be({s_v_end, s_v_id}) );
  s.push( to_be({pgt::complement(s_v_end), s_v_id}) );
  visited[s_v_idx] = 1;
  add_vertex_to_tree(s_v_end, s_v_idx);

  /* ---------- Main Loop ---------- */

  while (!s.empty()) {
    found_new_neighbour = false;
    pt::idx_t be_v_idx = s.top();

    p_idx = be_idx_to_ctr[be_v_idx];
    auto [syd, v_id] = to_bd(be_v_idx);
    pt::idx_t bd_v_idx = g.v_id_to_idx(v_id);

    const Vertex &v = g.get_vertex_by_id(v_id);
    const std::set<pt::idx_t> &neighbours = syd == pgt::v_end_e::l ? v.get_edges_l() : v.get_edges_r();

    // if no neighbours then it is a tip. Add a backedge to the root
    if (__builtin_expect((neighbours.empty() && !are_connected(p_idx, root_idx)), 0)) {
      t.add_be(p_idx, root_idx, pst::e_type_e::back_edge, pgt::color_e::gray);
      connect(p_idx, root_idx);
    }

    // stop processing edges on the first instance of finding a new neighbour
    for (auto e_idx : neighbours) {
      if ((found_new_neighbour = process_edge(bd_v_idx, syd, e_idx))) break;
    }

    if (!found_new_neighbour) s.pop();
  }

  return t;
}

} // namespace povu::bidirected
