#include "./bidirected.hpp"
#include "spanning_tree.hpp"
#include <cstddef>
#include <cstdint>
#include <iostream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>



namespace povu::bidirected {



/*
  Edge
  ------
 */
Edge::Edge(pt::idx_t v1_idx, pgt::v_end v1_end , pt::idx_t v2_idx, pgt::v_end v2_end)
  : v1_idx_{v1_idx}, v1_end{v1_end}, v2_idx_{v2_idx}, v2_end{v2_end} {}
pt::idx_t Edge::get_v1_idx() const { return this->v1_idx_; }
pgt::v_end Edge::get_v1_end() const { return this->v1_end; }
pt::idx_t Edge::get_v2_idx() const { return this->v2_idx_; }
pgt::v_end Edge::get_v2_end() const { return this->v2_end; }
pgt::side_n_id_t Edge::get_other_vertex(pt::idx_t v_idx) const {
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
  this->dummy_idx_ = constants::UNDEFINED_IDX;
}

/* getters */

pt::id_t VG::v_idx_to_id(pt::idx_t v_idx) const {
  return this->v_id_to_idx_.get_key(v_idx);
}

pt::idx_t VG::v_id_to_idx(pt::id_t v_id) const {
  return this->v_id_to_idx_.get_value(v_id);
}

pt::idx_t VG::size() const { return this->vertices.size(); }

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



void VG::add_tip(std::size_t v_id, pgt::v_end end) {
  this->tips_.insert( pgt::side_n_id_t{end, v_id} );
}

/* setters */

pt::idx_t VG::add_vertex(pt::id_t v_id, const std::string &label) {
  vertices.push_back(Vertex{v_id, label});
  this->v_id_to_idx_.insert(v_id, vertices.size() - 1);
  return vertices.size() - 1;
}

pt::idx_t VG::add_edge(pt::id_t v1_id, pgt::v_end v1_end, pt::id_t v2_id, pgt::v_end v2_end) {
  pt::idx_t v1_idx = this->v_id_to_idx_.get_value(v1_id);
  pt::idx_t v2_idx = this->v_id_to_idx_.get_value(v2_id);
  edges.push_back(Edge{v1_idx, v1_end, v2_idx, v2_end});
  pt::idx_t e_idx = edges.size() - 1;

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

  return e_idx;

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
  std::cout << "\t" << "vertex count: " << this->vtx_count() << std::endl;
  std::cout << "\t" << "edge count: " << this->edge_count() << std::endl;
  std::cout << "\t" << "Tip count " << this->tips().size() << std::endl;
  std::cout << "\t";
  pu::print_with_comma(std::cout, this->tips(), ',');
  std::cout << std::endl;
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

    os << "\t"
       << v_idx
       << std::format("[style=filled, fillcolor=lightblue, label=\"+ {} -\"];", v_id)
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

  bool found_unvisited { false };

  std::set<pt::idx_t> comp_vtxs; // current component vertices
  comp_vtxs.insert(start_vtx);

  /* ---------- Helper Functions ---------- */

  auto process_edge = [&](pt::idx_t v_idx, pt::idx_t e_idx) -> bool {
    const Edge &e = g.get_edge(e_idx);
    auto [_, adj_v_idx] = e.get_other_vertex(v_idx);
    if (adj_v_idx != v_idx && !visited.contains(adj_v_idx)) {
      s.push(adj_v_idx);
      visited.insert(adj_v_idx);
      comp_vtxs.insert(adj_v_idx);
      return true;
    }
    return false;
  };

  auto add_edges = [&](pt::idx_t v_idx, pt::idx_t e_idx, const Vertex &v) -> void {
    if (added_edges.contains(e_idx)) {
      return;
    }

    added_edges.insert(e_idx);
    const Edge &e = g.get_edge(e_idx);
    auto [side, adj_v_idx] = e.get_other_vertex(v_idx);
    curr_vg->add_edge(v.id(), pgt::v_end::r, g.v_idx_to_id(adj_v_idx), side);
  };

  /* ---------- Main Component Search Loop ---------- */

  while (!s.empty()) {
    found_unvisited = false;
    pt::idx_t v_idx = s.top();
    const Vertex& v = g.get_vertex_by_idx(v_idx);

    for (auto e_idx : v.get_edges_l()){
      found_unvisited = process_edge(v_idx, e_idx);
    }

    if (found_unvisited) {
      continue;
    }

    for (auto e_idx : v.get_edges_r()) {
      found_unvisited = process_edge(v_idx, e_idx);
    }

    if (!found_unvisited) {
      s.pop();
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
          add_edges(v_idx, e_idx, v);
        }

        for (auto e_idx : v.get_edges_r()) {
          add_edges(v_idx, e_idx, v);
        }
      }

      /* add tips */
      for (auto [side, v_id] : g.tips()) {
        if ( comp_vtxs.contains(g.v_id_to_idx(v_id))) {
          curr_vg->add_tip(v_id, side);
        }
      }

      // clear the set for the next component
      added_edges.clear();
      components.push_back(curr_vg);
      curr_vg = nullptr;

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

  /* assumes that the graph has a dummy vertex */
  if (g.get_dummy_idx() == constants::UNDEFINED_IDX) {
    std::cerr
      << "Graph has no dummy vertex computation of spanning tree expects dummy"
      << g.v_idx_to_id(0) << " "
      << std::endl;
    return pst::Tree { 0 }; // empty tree
  }

  // we bi-edge all vertices except for the single dummy start vertex
  pt::idx_t tree_size = ((g.vtx_count() - 1) * 2) + 1;
  pst::Tree t { tree_size };

  std::stack<side_n_id_t> s; // stores the vertex index not id

  // TODO: replace visited and v_idx_to_tree_idx with a single vector
  //       mapping side and id to a 2*idx + 1 value

  std::set<pgt::side_n_id_t> visited;
  std::map<side_n_id_t, pt::idx_t> v_idx_to_tree_idx;

  side_n_id_t start = { pgt::v_end::r, pc::DUMMY_VTX_ID };

  std::set<std::pair<id_t, id_t>> added_edges;

  s.push(start);
  visited.insert(start);
  // workaround for dummy vertex being bidirected
  visited.insert( {pgt::v_end::l, pc::DUMMY_VTX_ID} );

  bool found_unvisited { false };

  pt::idx_t counter{0};

  /* ---------- Helper Functions ---------- */

  auto mark_edge = [&](const side_n_id_t &a, const side_n_id_t &b) -> void {
    added_edges.insert({v_idx_to_tree_idx[a], v_idx_to_tree_idx[b]});
    added_edges.insert({v_idx_to_tree_idx[b], v_idx_to_tree_idx[a]});
  };

  auto is_marked = [&](const side_n_id_t a, const side_n_id_t &b) -> bool {
    return added_edges.contains( {v_idx_to_tree_idx[a], v_idx_to_tree_idx[b]} );
  };

  auto get_alt = [](pgt::v_end s) -> pgt::v_end {
    return s == pgt::v_end::l ? pgt::v_end::r : pgt::v_end::l;
  };

  auto end2typ = [](pt::id_t v_id, pst::v_end e) -> v_type_e {
    if (v_id == pc::DUMMY_VTX_ID) {
      return v_type_e::dummy;
    }
    return e == pgt::v_end::l ? v_type_e::l : v_type_e::r;
  };

  auto add_vertex_to_tree = [&](const side_n_id_t &to_add) -> void {
    if (!v_idx_to_tree_idx.contains(to_add)) {
      auto [end, v_id] = to_add;
      t.add_vertex({counter++, v_id, end2typ(v_id, end)});
      v_idx_to_tree_idx[to_add] = counter - 1;
    }
  };

  // returns true if a new vertex is discovered, false otherwise
  auto process_edge = [&](const side_n_id_t &curr, pt::idx_t e_idx) -> bool {
    auto [syd, v_id] = curr;
    pt::idx_t v_idx = g.v_id_to_idx(v_id);

    const Edge &e = g.get_edge(e_idx);
    auto [n_adj_side, n_v_idx] = e.get_other_vertex(v_idx);

    pt::id_t n_v_id = g.v_idx_to_id(n_v_idx);
    side_n_id_t neighbour = {n_adj_side, n_v_id};

    v_end_t n_alt_side = get_alt(n_adj_side);
    side_n_id_t n_alt = {n_alt_side, n_v_id};

    if ( !visited.contains(neighbour) ) {

      // when we discover a new vertex we create two vertices in the
      // spanning tree and join them with a black edge. We also add a a gray
      // tree edge between the parent and the child

      add_vertex_to_tree(curr);
      add_vertex_to_tree(neighbour); // discovered vertex
      add_vertex_to_tree(n_alt); // alt side of the discovered vertex

      t.add_tree_edge(v_idx_to_tree_idx[curr], v_idx_to_tree_idx[neighbour], color::gray);
      t.add_tree_edge(v_idx_to_tree_idx[neighbour], v_idx_to_tree_idx[n_alt], color::black);

      // DFS treversal state update

      mark_edge(curr, neighbour);
      mark_edge(neighbour, n_alt);

      s.push(neighbour);
      s.push(n_alt);

      visited.insert(neighbour);
      visited.insert(n_alt);

      return true; // new vertex discovered
    }
    else if (!is_marked(curr, neighbour)) {
      // add a backedge if:
      //  - not a parent child relationship
      //  - a backedge does not already exist
      t.add_be(v_idx_to_tree_idx[curr], v_idx_to_tree_idx[neighbour],
               pst::EdgeType::back_edge, color::gray);
      mark_edge(curr, neighbour);
    }

    return false; // no new vertex discovered
  };

  /* ---------- Main Loop ---------- */
  while (!s.empty()) {
    found_unvisited = false;
    side_n_id_t curr = s.top();
    auto [syd, v_id] = curr;

    const Vertex &v = g.get_vertex_by_id(v_id);
    const std::set<pt::idx_t>& neighbours = syd == pgt::v_end::l ? v.get_edges_l() : v.get_edges_r();

    for (auto e_idx : neighbours) {
      if (process_edge(curr, e_idx)) {
        found_unvisited = true;
        break;
      }
    }

    if (!found_unvisited) {
      s.pop();
    }
  }

  return t;
}

} // namespace povu::bidirected
