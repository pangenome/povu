#include "./bidirected.hpp"

namespace povu::bidirected {

// ============================================================
//      RefInfo
// ============================================================

// --------------
// constructor(s)
// --------------
RefInfo::RefInfo(pt::id_t ref_id, pgt::or_e strand, pt::idx_t locus)
    : ref_id_(ref_id), strand_(strand), locus_(locus) {}

// ---------
// getter(s)
// ---------
pt::id_t RefInfo::get_ref_id() const { return ref_id_; }
pgt::or_e RefInfo::get_strand() const { return strand_; }
pt::idx_t RefInfo::get_locus() const { return locus_; }

// ============================================================
//      Edge
// ============================================================


Edge::Edge(pt::idx_t v1_idx, pgt::v_end_e v1_end , pt::idx_t v2_idx, pgt::v_end_e v2_end)
  : v1_idx_{v1_idx}, v1_end_{v1_end}, v2_idx_{v2_idx}, v2_end_{v2_end} {}
pt::idx_t Edge::get_v1_idx() const { return this->v1_idx_; }
pt::idx_t &Edge::get_v1_idx_mut() { return this->v1_idx_; }
pgt::v_end_e Edge::get_v1_end() const { return this->v1_end_; }
pt::idx_t Edge::get_v2_idx() const { return this->v2_idx_; }
pt::idx_t &Edge::get_v2_idx_mut() { return this->v2_idx_; }
pgt::v_end_e Edge::get_v2_end() const { return this->v2_end_; }
pgt::side_n_idx_t Edge::get_other_vtx(pt::idx_t v_idx) const {
  return (get_v1_idx() == v_idx) ? pgt::side_n_id_t{v2_end_, v2_idx_}
                                 : pgt::side_n_id_t{v1_end_, v1_idx_};
}

pgt::side_n_idx_t Edge::get_other_vtx(pt::idx_t v_idx, pgt::v_end_e ve) const {
  const pt::idx_t v1 = get_v1_idx();
  const pt::idx_t v2 = get_v2_idx();

  if (v1 == v2) { // Handle self-loop case
    return {pgt::complement(ve), v1};
  }

  // Return the opposite vertex
  return (v1 == v_idx) ? pgt::side_n_id_t{v2_end_, v2_idx_} : pgt::side_n_id_t{v1_end_, v1_idx_};
}


/*
  Vertex
  ------
 */

// --------------
// constructor(s)
// --------------

Vertex::Vertex(pt::id_t v_id, const std::string& label) : v_id_{v_id}, label_(label) {}


// ---------
// getter(s)
// ---------

pt::id_t Vertex::id() const { return v_id_; }
const std::string &Vertex::get_label() const { return this->label_; }
std::string Vertex::get_rc_label() const {
  return pu::reverse_complement(this->label_);
}
const std::set<pt::idx_t>& Vertex::get_edges_l() const { return e_l; }
const std::set<pt::idx_t>& Vertex::get_edges_r() const { return e_r; }
const std::vector<RefInfo>& Vertex::get_refs() const { return refs_; }

// ---------
// setter(s)
// ---------

void Vertex::add_edge_l(pt::idx_t e_idx) { e_l.insert(e_idx); }
void Vertex::add_edge_r(pt::idx_t e_idx) { e_r.insert(e_idx); }
void Vertex::add_ref(pt::idx_t ref_id, pgt::or_e strand, pt::idx_t locus) {
  this->refs_.push_back(RefInfo(ref_id, strand, locus));
}


/*
  Graph
  -----
 */


// --------------
// constructor(s)
// --------------

VariationGraph::VariationGraph(pt::idx_t v_count, pt::idx_t e_count, bool inc_refs) {
  this->vertices.reserve(v_count);
  this->edges.reserve(e_count);

  this->has_refs_ = false;
  if (inc_refs) {
    this->has_refs_ = true;
    this->refs_ = Refs();
  }
}

// ---------
// getter(s)
// ---------

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
const Edge& VG::get_edge(pt::idx_t e_idx) const { return edges[e_idx]; }
Edge& VG::get_edge_mut(pt::idx_t e_idx) { return edges[e_idx]; }
const Vertex& VG::get_vertex_by_idx(pt::idx_t v_idx) const { return vertices[v_idx]; }
const Vertex& VG::get_vertex_by_id(pt::id_t v_id) const {
  return vertices[this->v_id_to_idx_.get_value(v_id)];
}
Vertex& VG::get_vertex_mut_by_id(pt::id_t v_id) {
  return vertices[this->v_id_to_idx_.get_value(v_id)];
}
const std::string &VG::get_ref_label(pt::id_t ref_id) const {
  return this->refs_.get_ref_label(ref_id);
}

const pgt::Ref &VG::get_ref_by_id(pt::id_t ref_id) const {
  return this->refs_.get_ref(ref_id);
}

pgt::Ref &VG::get_ref_by_id_mut(pt::id_t ref_id) {
  return this->refs_.get_ref_mut(ref_id);
}

pt::id_t VG::get_ref_id(const std::string &ref_label) const {
  return this->refs_.get_ref_id(ref_label);
}

// ---------
// setter(s)
// ---------

void VG::add_tip(pt::id_t v_id, pgt::v_end_e end) {
  this->tips_.insert(pgt::side_n_id_t{end, v_id});
}

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

pt::id_t VG::add_ref(const std::string &label, char delim) {
  std::string fn_name{pv_cmp::format("[{}::VG::{}]", MODULE, __func__)};
  pt::id_t ref_id = this->refs_.add_ref(label, delim);

  return ref_id;
}

const std::set<pt::id_t> &VG::get_shared_samples(pt::id_t ref_id) const {
  return this->refs_.get_shared_samples(ref_id);
}

pt::id_t VG::ref_id_count() const {
  return this->refs_.ref_id_count();
}


void VG::shrink_to_fit() {
  this->vertices.shrink_to_fit();
  this->edges.shrink_to_fit();
}

void VG::summary(bool print_tips) const {
  std::cout << "Bidirected Graph: " << std::endl;
  std::cout << "\t" << "vertex count: " << this->vtx_count() << std::endl;
  std::cout << "\t" << "edge count: " << this->edge_count() << std::endl;
  std::cout << "\t" << "Tip count " << this->tips().size() << std::endl;
  if (print_tips) {
    std::cerr << "\t" << "Tips: ";
    std::cout << "\t";
    pu::print_with_comma(std::cout, this->tips(), ',');
    std::cout << std::endl;
  }
}

void VG::print_dot(std::ostream& os) const {

  const std::string header = R"(
graph G {
    graph [rankdir=LR];
    node [shape=cds, style=filled, fillcolor=lightblue, fontsize="10pt"];
)";

  /* helper fns */
  // map v end left and right to dot west and east for rectangular vertices
  auto v_end_to_dot = [](pgt::v_end_e e) -> std::string {
    return e == pgt::v_end_e::l ? "w" : "e";
  };

  /* header */
  os << header;

  /* vertices */
  for (size_t v_idx {}; v_idx < this->vtx_count(); ++v_idx) {
    const Vertex& v = this->get_vertex_by_idx(v_idx);
    std::string v_id = v.id() == constants::UNDEFINED_ID ? "d" : std::to_string(v.id());

    os << pv_cmp::format("\t{}[label=\"+ {} - \\n ({})\"];\n", v_idx, v_id, v_idx);
  }

  /* edges */
  for (const Edge& e: this->edges) {
    pt::idx_t v1_idx = e.get_v1_idx();
    std::string v1_e = v_end_to_dot(e.get_v1_end());
    pt::idx_t v2_idx = e.get_v2_idx();
    std::string v2_e = v_end_to_dot(e.get_v2_end());

    os << pv_cmp::format("\t{}:{}--{}:{}[color=gray];\n", v1_idx, v1_e, v2_idx, v2_e);
  }

  /* footer */
  os << "}" << std::endl;
}

// TODO: make this static factory fn
// does not handle refs, should it?
std::vector<VG *> componetize(const povu::bidirected::VG &g) {

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

    if (pv_cmp::contains(visited, adj_v_idx)) { // also handles self loops
      return;
    }

    s.push(adj_v_idx);
    visited.insert(adj_v_idx);
    comp_vtxs.insert(adj_v_idx);
    return;
  };

  auto add_edges = [&](const Vertex &v, pgt::v_end_e ve, pt::idx_t v_idx, pt::idx_t e_idx) -> void {
   if (pv_cmp::contains(added_edges, e_idx)) { // don't duplicate edges
       return;
   }

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
      curr_vg = new VG(comp_vtxs.size(), added_edges.size(), false);

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
        if (pv_cmp::contains(comp_vtxs, g.v_id_to_idx(v_id))) {
          curr_vg->add_tip(v_id, side);
        }
      }

      // clear the set for the next component
      components.push_back(curr_vg);
      curr_vg = nullptr;
      added_edges.clear();

      /* find the next unvisited vertex */
      for (std::size_t v_idx{}; v_idx < g.vtx_count(); ++v_idx) {
        if (!pv_cmp::contains(visited, v_idx)) { // if not visited
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

} // namespace povu::bidirected
