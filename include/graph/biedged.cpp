#include "./biedged.hpp"
#include <cstddef>

namespace biedged {

/*
  ====
  Edge
  ====
*/
Edge::Edge(std::size_t v1_idx, v_type v1_type, std::size_t v2_idx, v_type v2_type, pgt::color_e c)
  : v1_idx(v1_idx), v1_type(v1_type), v2_idx(v2_idx), v2_type(v2_type), c(c)
{
  std::string fn_name = "[povu::graph::biedged::Edge]";

  if (v1_idx == v2_idx) {
    throw std::invalid_argument(std::format("{} Self-loops are not allowed Biedged vertex indexes {} {}", fn_name, v1_idx, v2_idx));
  }
}

std::size_t Edge::get_v1_idx() const {
    return this->v1_idx;
}

std::size_t Edge::get_v2_idx() const {
    return this->v2_idx;
}

color Edge::get_color() const {
  return this->c;
}

std::size_t Edge::get_eq_class() const {
  return this->eq_class;
}

std::size_t Edge::get_other_vertex(std::size_t vertex_index) const {
  std::string fn_name = std::format("[povu::biedged::{}]", __func__);
  if (vertex_index == this->v1_idx) {
    return this->v2_idx;
  }
  else if (vertex_index == this->v2_idx) {
    return this->v1_idx;
  }
  else {
    throw std::invalid_argument(std::format("{} Vertex {} is not part of edge",
                                            fn_name, vertex_index));
  }
}



void Edge::set_eq_class(std::size_t eq_class) {
  this->eq_class = eq_class;
}

void Edge::set_v1_idx(std::size_t i) {
    this->v1_idx = i;
}

void Edge::set_v2_idx(std::size_t i) {
    this->v2_idx = i;
}

std::ostream& operator<<(std::ostream& os, const Edge& e) {
  os << "Edge (" << e.v1_idx << ", " << e.v1_type << e.v2_idx << ", " << e.v2_type << e.c << ")";
  return os;
}

/*
  ======
  Vertex
  ======
*/

/*
  Constructor(s)
*/

Vertex::Vertex() :
  black_edge(pc::UNDEFINED_SIZE_T),
  grey_edges(std::set<std::size_t>()),
  type(v_type::l),
  vertex_idx(0)
{}

Vertex::Vertex(std::size_t id, std::size_t vertex_idx, v_type vertex_type) :
  black_edge(pc::UNDEFINED_SIZE_T),
  grey_edges(std::set<std::size_t>()),
  type(vertex_type),
  id_(id),
  vertex_idx(vertex_idx)
{}


std::size_t Vertex::id() const {
    return this->id_;
}

std::set<std::size_t> const& Vertex::get_grey_edges() const {
    return this->grey_edges;
}

std::size_t Vertex::get_black_edge() const {
    return this->black_edge;
}

std::vector<std::size_t> Vertex::get_neighbours() const {
  std::vector<std::size_t> neighbours;
  if (this->get_black_edge() != pc::UNDEFINED_SIZE_T) {
    neighbours.push_back(this->get_black_edge());
  }
  // add all grey edges to the neighbours vector
  for (auto const& edge : this->get_grey_edges()) {
    neighbours.push_back(edge);
  }

  return neighbours;
}

std::size_t Vertex::get_vertex_idx() const {
  return this->vertex_idx;
}

void Vertex::add_edge(std::size_t edge_idx, color c) {
  if (c == color::black) {
    black_edge = edge_idx;
  } else {
    grey_edges.insert(edge_idx);
  }
}

v_type Vertex::get_type() const {
    return type;
}

void Vertex::set_vertex_idx(std::size_t i) {
    this->vertex_idx = i;
}

/*
  =========================================
  BVariationGraph (Biedged Variation Graph)
  =========================================
*/

/*
  Constructor(s)
*/
void validate_biedging(BVariationGraph const& bg) {
  for (std::size_t i{} ; i < bg.size(); ++i) {
    if (bg.get_neighbours(i).size() < 2) {
      std::cout << "failed check i " << i << std::endl;
    }

    if (i>0 && i < bg.size() - 1) {
      if (bg.get_vertex(i).get_black_edge() == pc::UNDEFINED_SIZE_T) {
        std::cout << "no black edge failed check i " << i << std::endl;
      }
    }
  }
}



BVG::BVariationGraph(const pbd::VG &g) {
  std::string fn_name = std::format("[povu::BVariationGraph::{}]", __func__);

  std::size_t biedged_size = g.size() * 2 + 1;

  this->vertices = std::vector<biedged::Vertex>();
  this->vertices.reserve(biedged_size);

  /* lambdas to map a vertex at index i in the bidirected graph to an index in the biedged graph */
  // 1 is the dummy vertex count
  // TODO: replace or offload to the already existing fns in utils namespace
  auto new_l = [&](std::size_t i) { return (2*(i+1)) - 2 + 1 ; };
  auto new_r = [&](std::size_t i) { return (2*(i+1)) - 1 + 1 ; };

  /* add a dummy vertex */
  // given GFA starts at 1 let the dummy vertex id be 0
  const std::size_t d_vtx_idx = 0;
  const std::size_t d_vtx_id = 0;
  this->vertices.push_back(Vertex(d_vtx_id, d_vtx_idx, v_type::dummy));

  /* add vertices and black edges */
  std::size_t l_idx, r_idx, v_id;
  for (size_t v_idx {}; v_idx < g.size(); ++v_idx) {
    l_idx = new_l(v_idx);
    r_idx = new_r(v_idx);
    v_id = g.v_idx_to_id(v_idx);

    this->vertices.push_back(Vertex(v_id, l_idx, v_type::l));
    this->vertices.push_back(Vertex(v_id, r_idx, v_type::r));
    this->add_edge(l_idx, v_type::l, r_idx, v_type::r, color::black);
  }


  /* add grey edges */
  std::size_t v1_idx, v2_idx;
  for(std::size_t e_idx {}; e_idx < g.edge_count(); e_idx++) {
    const pbd::Edge &e = g.get_edge(e_idx);
    v1_idx = e.get_v1_idx();
    v2_idx = e.get_v2_idx();
    auto [l_v1_idx, r_v1_idx] = pu::frm_bidirected_idx(v1_idx);
    auto [l_v2_idx, r_v2_idx] = pu::frm_bidirected_idx(v2_idx);

    if (e.get_v1_end() == pgt::v_end::l && e.get_v2_end() == pgt::v_end::l) { // l l
      this->add_edge(l_v1_idx, v_type::l, l_v2_idx, v_type::l, pgt::color_e::gray);
    }
    else if (e.get_v1_end() == pgt::v_end::r && e.get_v2_end() == pgt::v_end::l) { // r l
      this->add_edge(r_v1_idx, v_type::r, l_v2_idx, v_type::l, pgt::color_e::gray);
    }
    else if (e.get_v1_end() == pgt::v_end::l && e.get_v2_end() == pgt::v_end::r) { // l r
      this->add_edge(l_v1_idx, v_type::l, r_v2_idx, v_type::r, pgt::color_e::gray);
    }
    else { // r r
      this->add_edge(r_v1_idx, v_type::r, r_v2_idx, v_type::r, pgt::color_e::gray);
    }
  }

  /* connect dummy vertex to tips  */
  Vertex be_v;
  for (std::size_t be_v_idx{}; be_v_idx < biedged_size; be_v_idx++ ) {
    be_v = this->get_vertex_mut(be_v_idx);
    if (be_v.get_type() == v_type::dummy) { continue; }
    if (be_v.get_grey_edges().empty()) {
      this->add_edge(d_vtx_id, v_type::dummy, be_v_idx, be_v.get_type(), color::gray);
    }
  }

#ifdef DEBUG
  validate_biedging(*this);
#endif
}


/*
  Getters
*/
std::vector<std::pair<color, std::size_t>> BVariationGraph::get_neighbours(std::size_t vertex_idx) const {
  std::vector<std::pair<color, std::size_t>> neighbours;
  Vertex const& v = this->get_vertex(vertex_idx);


  std::size_t black_e_idx = v.get_black_edge();
  if (black_e_idx != pc::UNDEFINED_SIZE_T) {
  neighbours.push_back(std::make_pair(
      color::black,
      this->get_edge(black_e_idx).get_other_vertex(vertex_idx)));
  }


  for (auto e_idx : v.get_grey_edges()) {
    const Edge& e = this->get_edge(e_idx);
    std::size_t other_vertex = e.get_other_vertex(vertex_idx);
    std::pair<color, std::size_t> p = std::make_pair(color::gray, other_vertex);
    neighbours.push_back(p);
  }


  return neighbours;
}

const biedged::Vertex& BVariationGraph::get_vertex(std::size_t i) const {
  return vertices.at(i);
}

const std::vector<size_t>& BVariationGraph::get_dummy_vertices() const {
  return this->dummy_vertices_;
}

std::size_t BVariationGraph::size() const {
  return this->vertices.size();
}

std::size_t BVariationGraph::num_edges() const {
  return this->edges.size();
}


Vertex& BVariationGraph::get_vertex_mut(std::size_t i) {
  return vertices.at(i);
}

const std::vector<Edge>& BVariationGraph::get_all_edges() const {
  return this->edges;
}

Edge& BVariationGraph::get_edge_mut(std::size_t e_idx) {
  return this->edges.at(e_idx);
}

const Edge& BVariationGraph::get_edge(std::size_t edge_idx) const {
  return this->edges.at(edge_idx);
}


/*
  Setters
*/

std::size_t BVariationGraph::add_edge(std::size_t v1_idx, v_type v1_type,
                                      std::size_t v2_idx, v_type v2_type,
                                      color c) {
  this->edges.push_back(Edge(v1_idx, v1_type, v2_idx, v2_type, c));
  this->get_vertex_mut(v1_idx).add_edge(this->edges.size() - 1, c);
  this->get_vertex_mut(v2_idx).add_edge(this->edges.size() - 1, c);

  return this->edges.size() - 1;
}


void BVariationGraph::add_vertex(std::size_t v_id, v_type vertex_type) {
  vertices.push_back(Vertex(v_id, vertices.size(), vertex_type));
}


void BVariationGraph::update_eq_classes(pst::Tree& st) {
  const std::map<std::size_t, std::pair<pst::EdgeType, std::size_t>>& m = st.get_g_edge_idx_map();

  for (auto& [e_idx, v] : m) {
    auto [e_type, t_e_idx] = v;

    std::size_t eq_class = (e_type == pst::EdgeType::tree_edge)
      ? st.get_tree_edge(v.second).get_class()
      : st.get_backedge(v.second).get_class();

    //auto e = this->get_edge_mut(e_idx);

    this->get_edge_mut(e_idx).set_eq_class(eq_class);
  }
}

/*
  Misc
*/

void BVariationGraph::print_dot() const {
  std::cout <<
    "graph G {\n" <<
    "\trankdir=LR;\n" <<
    "\tnode [shape=circle];\n";

  for (std::size_t v_idx{}; v_idx < this->vertices.size(); ++v_idx) {
    std::cout <<
      std::format("\t{} [label=\"{}{}\"]\n",
                  v_idx,
                  this->get_vertex(v_idx).id(),
                  (this->get_vertex(v_idx).get_type() == v_type::l ? "+" : "-") );
  }

  //
  for (std::size_t e_idx{}; e_idx < this->edges.size() ; ++e_idx) {
    auto e = this->edges.at(e_idx);
    if (e.get_color() == color::black) {
      std::cout <<
        std::format("\t{} -- {} [color=\"black\"; label=\"{}\"];\n",
                    e.get_v1_idx(), e.get_v2_idx(), (e.get_eq_class() == pc::UNDEFINED_SIZE_T ? "" : std::to_string(e.get_eq_class())) );
    }
    else {
      std::cout <<
        std::format("\t{} -- {} [color=\"gray\" label=\"{}\"];\n", e.get_v1_idx(), e.get_v2_idx(),
                    (e.get_eq_class() == pc::UNDEFINED_SIZE_T ? "" : std::to_string(e.get_eq_class())));
    }
  }

  std::cout << "}" << std::endl;
}

// TODO: make ref params const
void validate_spanning_tree(BVariationGraph const &g, pst::Tree &t,
                            //const pu::TwoWayMap<std::size_t, std::size_t>& vtx_to_dfs_num
                            std::map<std::size_t, std::size_t>& dfs_num_to_vtx) {
  if (t.size() != g.size()) {
    std::cerr << "size mismatch " << t.size() << " " << g.size() << "\n";
  }

  // compare edge count
  if (t.tree_edge_count() + t.back_edge_count() != g.num_edges()) {
    std::cerr << std::format("edge count mismatch {} {} {}\n",
                             t.tree_edge_count(), t.back_edge_count(), g.num_edges());
  }

  for (std::size_t v{}; v < t.size() ; v++) {
    std::set<std::size_t>const& children = t.get_children(v);
    std::set<std::size_t>const& obe = t.get_obe(v);
    std::set<std::size_t> const &ibe = t.get_ibe(v);

    //dfs_num_to_vtx[v];

    if (t.is_root(v)) {
      if (obe.size() + ibe.size() + children.size() != g.get_neighbours(    dfs_num_to_vtx[v]).size()) {
        std::cerr << std::format("neighbour mismatch (root) {} {} {} {} {}\n",
                                 v,
                                 obe.size(),
                                 ibe.size(),
                                 children.size(),
                                 g.get_neighbours(v).size());
      }
    }
    else {
      if (obe.size() + ibe.size() + children.size() + 1 != g.get_neighbours(    dfs_num_to_vtx[v]).size()) {
        std::cerr << std::format("neighbour mismatch {} ({} {}) {} {} {} {}\n",
                                 v,
                                 t.get_vertex(v).g_v_id(),
                                 g.get_vertex(dfs_num_to_vtx[v]).id(),
                                 obe.size(),
                                 ibe.size(),
                                 children.size(), g.get_neighbours(v).size());
      }
    }

    for (auto c : t.get_children(v)) {
      if (t.get_vertex(v).dfs_num() >= t.get_vertex(c).dfs_num()) {
        g.get_vertex( dfs_num_to_vtx[v]).id();
        std::cerr << std::format("te weird {} {} {} {}\n",
                                 v,
                                 t.get_vertex(v).dfs_num(),
                                 c,
                                 t.get_vertex(c).dfs_num());
      }
    }

    for (std::size_t tgt : t.get_obe(v)) {
      if (t.get_vertex(v).dfs_num() <= t.get_vertex(tgt).dfs_num()) {
        std::cerr << std::format("obe weird {} {} {} {}\n",
                                 v,
                                 t.get_vertex(v).dfs_num(),
                                 tgt,
                                 t.get_vertex(tgt).dfs_num());
      }
    }

    for (std::size_t src : t.get_ibe(v)) {
      if (t.get_vertex(v).dfs_num() >= t.get_vertex(src).dfs_num()) {
        std::cerr << std::format("ibe weird {} {} {} {}\n",
                                 v,
                                 t.get_vertex(v).dfs_num(),
                                 src,
                                 t.get_vertex(src).dfs_num());
      }
    }
  }
}


pst::Tree BVariationGraph::compute_spanning_tree() const {
  pst::Tree t = pst::Tree(this->size());

  std::stack<std::tuple<std::size_t, std::size_t, pgt::color>> s; // parent indedx in the biedged graph, vertex idx and color

  std::vector<bool> visited(this->size(), false); // visited vertices

  std::size_t v_idx {}; // set start node to 0
  s.push({pc::INVALID_ID, v_idx, color::gray});
  visited[v_idx] = true;

  std::size_t counter {}; // dfs pre visit counter

  std::map<std::size_t, std::size_t> vtx_to_dfs_num;

  std::vector<bool> in_tree(this->size(), false); // graph vertices in the tree
  std::set<pt::unordered_pair<std::size_t>> back_edges; // avoids duplicate back edges

  auto is_parent = [&](std::size_t p, std::size_t c) -> bool {
    return t.get_vertex(vtx_to_dfs_num[c]).get_parent_idx() == vtx_to_dfs_num[p];
  };

  auto is_parent_child = [&](std::size_t n, std::size_t v_idx) -> bool {
    return is_parent(n, v_idx) || is_parent(v_idx, n);
  };

  auto is_self_loop = [&](std::size_t n, std::size_t v_idx) -> bool {
    if (!is_parent_child(n, v_idx)) {
      return false;
    }

    if (this->get_vertex(v_idx).get_type() == v_type::dummy || this->get_vertex(n).get_type() == v_type::dummy) {
      return false;
    }

    std::size_t v_handle = this->get_vertex(v_idx).id();
    std::size_t n_handle = this->get_vertex(n).id();

    if (n_handle != v_handle) {
      return false;
    }


    for (auto e_idx : this->get_vertex(v_idx).get_grey_edges()) {
      if (this->get_edge(e_idx).get_other_vertex(v_idx) == n) {
        return true;
      }
    }

    return false;
  };


  while (!s.empty()) {
    auto [p_idx, v_idx, c] = s.top();

    biedged::Vertex const& v = this->get_vertex(v_idx);

    if (!in_tree[v_idx]) {
      t.add_vertex(pst::Vertex{counter, v.id(), v.get_type()});
      in_tree[v_idx] = true;
      vtx_to_dfs_num[v_idx] = counter;

      if (p_idx != pc::INVALID_ID) {
        t.add_tree_edge(vtx_to_dfs_num[p_idx], counter, pc::UNDEFINED_SIZE_T, c);
      }

      counter++;
    }

    bool explored {true};
    for (auto [c, n] : this->get_neighbours(v_idx)) {
      if (!visited[n]) {
        s.push({v_idx, n, c});
        visited[n] = true;
        explored = false;
        break;
      }
      else if ((is_self_loop(n, v_idx) || !is_parent_child(n, v_idx)) && !back_edges.count({n, v_idx}))  {
        t.add_be(vtx_to_dfs_num[v_idx], vtx_to_dfs_num[n], pc::UNDEFINED_SIZE_T, pst::EdgeType::back_edge, c);
        back_edges.insert({v_idx, n});
      }
    }


    if (explored) { s.pop(); }
  }

  //#ifdef DEBUG   // check the correctness of the tree
  // validate_spanning_tree(*this, t, dfs_num_to_vtx);
  //#endif

  return t;
}

}; // namespace biedged
