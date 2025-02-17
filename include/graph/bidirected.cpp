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
std::string Vertex::get_rc_label() const {
  return pu::reverse_complement(this->label_);
}
const std::set<pt::idx_t>& Vertex::get_edges_l() const { return e_l; }
const std::set<pt::idx_t>& Vertex::get_edges_r() const { return e_r; }
const std::vector<PathInfo>& Vertex::get_refs() const { return refs_; }
/* setters */
void Vertex::add_edge_l(pt::idx_t e_idx) { e_l.insert(e_idx); }
void Vertex::add_edge_r(pt::idx_t e_idx) { e_r.insert(e_idx); }
void Vertex::add_ref(pt::idx_t path_id, pgt::or_e strand, pt::idx_t step_index) {
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

    os << std::format("\t{}[label=\"+ {} - \\n ({})\"];\n", v_idx, v_id, v_idx);
  }

  /* edges */
  for (const Edge& e: this->edges) {
    pt::idx_t v1_idx = e.get_v1_idx();
    std::string v1_e = v_end_to_dot(e.get_v1_end());
    pt::idx_t v2_idx = e.get_v2_idx();
    std::string v2_e = v_end_to_dot(e.get_v2_end());

    os << std::format("\t{}:{}--{}:{}[color=gray];\n", v1_idx, v1_e, v2_idx, v2_e);
  }

  /* footer */
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
      t.add_be(p_idx, be_idx_to_ctr[o_be_idx], pst::be_type_e::back_edge, pgt::color_e::gray);
      connect(p_idx, be_idx_to_ctr[o_be_idx]);
    }
    else if (__builtin_expect((bd_v_idx == ov_idx && !self_loops.contains(bd_v_idx)), 0)) {
      // add a self loop backedge, a parent-child relationship
      t.add_be(p_idx, be_idx_to_ctr[o_be_idx] , pst::be_type_e::back_edge, pgt::color_e::gray);
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
      t.add_be(p_idx, root_idx, pst::be_type_e::back_edge, pgt::color_e::gray);
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


void populate_walks(const VG &g, pvt::RoV &r, pt::idx_t max_steps) {
  const std::string fn_name = std::format("[povu::bidirected::{}]", __func__);

  typedef id_or_t idx_or_t; // specifically for idx instead of id
  typedef idx_or_t step;
  enum class dir_e { in, out }; // direction

  std::queue<idx_or_t> q;

  // a map to keep track of the vertices whose incoming neighbours paths we have
  // extended so far key is the vertex and value is the set of vertices whose
  // paths we have extended
  std::map<idx_or_t, std::set<idx_or_t>> seen;

  // a set to keep track of the vertices we've seen
  std::set<idx_or_t> explored;

  bool all_incoming_explored{true};

  std::size_t counter {}; // a counter to keep track of the number of iterations
  // allows us to short circuit the traversal if counter > max_steps

  // a map of an idx and orientation to the incoming walk from the entry to the
  // idx and side that is incoming for the orientation
  std::map<idx_or_t, std::vector<pvt::Walk>> in_walks;

  // Returns the edge set for a vertex based on orientation and direction.
  auto get_edges = [](const Vertex &v, pgt::or_e o, dir_e d) -> const std::set<pt::idx_t> & {
    // For "in" direction, use left edges if forward; right otherwise.
    // For "out" direction, swap the logic.
    return d == dir_e::in
               ? (o == pgt::or_e::forward ? v.get_edges_l() : v.get_edges_r())
               : (o == pgt::or_e::forward ? v.get_edges_r() : v.get_edges_l());
  };

  auto or_to_v_end = [](pgt::or_e o, dir_e d) -> pgt::v_end_e {
    return d == dir_e::in
      ? (o == pgt::or_e::forward ? pgt::v_end_e::l : pgt::v_end_e::r)
      : (o == pgt::or_e::forward ? pgt::v_end_e::r : pgt::v_end_e::l);
  };

  // Given a vertex end side, compute the alternate orientation.
  auto get_alt_or = [](pgt::v_end_e side, dir_e d) -> pgt::or_e {
    return d == dir_e::in
      ? (side == pgt::v_end_e::r ? pgt::or_e::forward : pgt::or_e::reverse)
      : (side == pgt::v_end_e::l ? pgt::or_e::forward : pgt::or_e::reverse);
  };

  auto get_neighbours = [&](idx_or_t idx_n_o, dir_e d) -> std::set<idx_or_t> {
    auto [v_idx, o] = idx_n_o;
    std::set<idx_or_t> neighbours;
    pgt::v_end_e ve = or_to_v_end(o, d);
    const Vertex &v = g.get_vertex_by_idx(v_idx);

    //pgt::or_e alt_o;
    for (const auto &e_idx : get_edges(v, o, d)) {
      const Edge &e = g.get_edge(e_idx);
      auto [side, alt_idx] = e.get_other_vtx(v_idx, ve);
      neighbours.insert({alt_idx, get_alt_or(side, d)});
    }

    return neighbours;
  };

  auto append_q = [&](idx_or_t current) {
    for (const idx_or_t &out_n : get_neighbours(current, dir_e::out)) {
      if (!explored.contains(current) || !explored.contains(out_n)) {
        q.push(out_n);
      }
    }
  };

  //const pgt::id_or_t &entry = r.get_entry();
  //const pgt::id_or_t &exit = r.get_exit();

  auto [start_id, start_o] = r.get_entry();
  auto [stop_id, stop_o] = r.get_exit();

  pt::idx_t start_idx = g.v_id_to_idx(start_id);
  pt::idx_t stop_idx = g.v_id_to_idx(stop_id);

  idx_or_t s = {start_idx, start_o};
  idx_or_t t = {stop_idx, stop_o};

  /* initialise the traversal */
  q.push(s);
  in_walks[{start_id, start_o}].emplace_back( pvt::Walk{start_id, start_o} );
  append_q({g.v_idx_to_id(start_id), start_o});

  while (!q.empty()) {
    if (counter++ > max_steps) {
      std::cerr << fn_name << " max_steps reached for flubble "
                << r.get_entry() << " ~> "
                << r.get_exit() << std::endl;
      break;
    }

    idx_or_t current = q.front();
    auto [c_v_idx, c_o] = current;
    pgt::id_or_t c_id_or { g.v_idx_to_id(c_v_idx), c_o };
    q.pop();

    all_incoming_explored = true;


    for (const idx_or_t &n : get_neighbours(current, dir_e::in)) {

      auto [n_idx, n_o] = n;
      // if we've added paths from this neighbour before

      // by default the start will be explored
      if (!explored.count(n) && n_idx != c_v_idx) {
        all_incoming_explored = false;
      }

      if (seen[current].contains(n) || !explored.contains(n)) {
        continue;
      }

      // append the current step to the incoming walks of n
      in_walks[current] = in_walks[n];
      for (pvt::Walk &w : in_walks[current]) {
        w.append_step(c_id_or);
      }

      seen[current].insert(n);
    }

    if (current != t) {
      append_q(current);
    }

    if (all_incoming_explored) {
      explored.insert(current);
    }

  }

  r.set_walks(std::move(in_walks[{stop_id, stop_o}]));
  return;
}

} // namespace povu::bidirected
