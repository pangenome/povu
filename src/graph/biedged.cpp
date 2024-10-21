#include <cstddef>
#include <format>
#include <iostream>
#include <set>
#include <stdexcept>
#include <string>
#include <vector>
#include <stack>

#include "./biedged.hpp"
#include "../common/utils.hpp"
#include "../common/types.hpp"


namespace biedged {
namespace pu = povu::utils;
namespace pt = povu::types;
namespace pgt = povu::graph_types;

/*
  ====
  Edge
  ====
*/
Edge::Edge(std::size_t v1, v_type v1_type,
           std::size_t v2, v_type v2_type,
           color c)
  : v1_idx(v1), v1_type(v1_type), v2_idx(v2), v2_type(v2_type), c(c)
{
  std::string fn_name = "[povu::graph::biedged::Edge]";

  if (v1 == v2) {
    throw std::invalid_argument(std::format("{} Self-loops are not allowed {} {}", fn_name, v1, v2));
  }

  if (c == color::black) {
    // throw an argument error if the edge is black and no label is provided
    throw std::invalid_argument(fn_name + " edges must have a label");
  }

  this->label = std::string();
}

Edge::Edge(std::size_t v1, v_type v1_type,
           std::size_t v2, v_type v2_type,
           color c, std::string label)
  : v1_idx(v1), v1_type(v1_type), v2_idx(v2), v2_type(v2_type), c(c), label(label)
{
  std::string fn_name = "[povu::graph::biedged::Edge]";

  if (v1 == v2) {
    throw std::invalid_argument(fn_name + " Self-loops are not allowed");
  }

  if (c == color::gray) {
    // throw an argument error if the edge is gray and a label is provided
    throw std::invalid_argument(fn_name + " Gray edges cannot have a label");
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
    throw std::invalid_argument(std::format("{} Vertex {} is not part of edge", fn_name, vertex_index));
  }
}

const std::string& Edge::get_label() const {
    return this->label;
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
  handle(std::string()),
  vertex_idx(0)
{}

Vertex::Vertex(const std::string& id, std::size_t vertex_idx, v_type vertex_type) :
  black_edge(pc::UNDEFINED_SIZE_T),
  grey_edges(std::set<std::size_t>()),
  type(vertex_type),
  handle(id),
  vertex_idx(vertex_idx)
{}


const std::string& Vertex::get_handle() const {
    return this->handle;
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

BVariationGraph::BVariationGraph(const povu::graph::Graph &g, bool add_dummy_vertices) {
  std::string fn_name = std::format("[povu::BVariationGraph::{}]", __func__);


  std::size_t biedged_size = g.size() * 2 + (add_dummy_vertices ? 2 : 0);

  this->vertices = std::vector<biedged::Vertex>();
  this->vertices.reserve(biedged_size);

  // this->vertices = std::vector<biedged::Vertex>(biedged_size, biedged::Vertex());

  auto new_l = [&](std::size_t i) { return (2*(i+1)) - 2 + (add_dummy_vertices ?  1 : 0 ); };
  auto new_r = [&](std::size_t i) { return (2*(i+1)) - 1 + (add_dummy_vertices ?  1 : 0 ); };

  auto duplicate_vertices = [this, g](std::size_t v_id, std::size_t i_l, std::size_t i_r) {

    // --------------------
    // add vertices + and -
    // --------------------

    std::string v_id_str = std::to_string(v_id); // TODO: stop passing a string for ID?
    {
      this->vertices.push_back(Vertex(v_id_str, i_l, v_type::l));
      this->vertices.push_back(Vertex(v_id_str, i_r, v_type::r));
    }

    // ----------------
    // add a black edge
    // ----------------
    {
      std::size_t e_idx =
        this->add_edge(i_l, v_type::l,
                       i_r, v_type::r,
                       color::black,
                       "");     // TODO: URGENT: remove label

      this->get_vertex_mut(i_l).add_edge(e_idx, color::black);
      this->get_vertex_mut(i_r).add_edge(e_idx, color::black);
    }
  };

  std::set<pt::unordered_pair<std::size_t>> added_edges;

  /*
    add gray edges
  */
  auto do_gray_edges = [&](std::size_t v_idx, povu::graph::Vertex const& current_vg_vertex, std::size_t i_l, std::size_t i_r) {
      // add gray edges incident with the 5' (left/+) vertex
      for (auto e_idx : current_vg_vertex.get_edges_l()) {
        const povu::graph::Edge& e = g.get_edge(e_idx);

        std::size_t new_v1{};
        v_type v1_type{};

        if (e.get_v1_idx() == v_idx) {
          auto [l, r] =  pu::frm_bidirected_idx(e.get_v2_idx());
          new_v1 = e.get_v2_end() == v_end::r ? r : l;
          v1_type = e.get_v2_end() == v_end::l ? v_type::r : v_type::l;
        }
        else {
          auto [l, r] =  pu::frm_bidirected_idx(e.get_v1_idx());
          new_v1 = e.get_v1_end() == v_end::l ? l : r;
          v1_type = e.get_v1_end() == v_end::l ? v_type::l : v_type::r;
        }

        if (added_edges.count(pt::unordered_pair(new_v1, i_l))) { continue; }

        added_edges.insert(pt::unordered_pair(new_v1, i_l));

        if (new_v1 == i_l) {
          std::cout << "error edge e: " << e.get_v1_idx() << std::endl;
          std::cout << "v_idx " << v_idx << " " << e.get_v1_idx() << " " << e.get_v2_end() << std::endl;
          std::cout << "n " << new_v1 << "r " << new_r(e.get_v2_idx()) << " l " << new_l(e.get_v2_idx()) << std::endl;
          //continue;
        }

        std::size_t gray_e_idx = this->add_edge(new_v1, v1_type, i_l, v_type::l,color::gray);

        //this->get_vertex_mut(new_v1).add_edge(gray_e_idx, color::gray);
        //this->get_vertex_mut(i_l).add_edge(gray_e_idx, color::gray);
      }

      // add gray edges incident with the 3' (right/-) vertex
      for (auto e_idx : current_vg_vertex.get_edges_r()) {

        const povu::graph::Edge& e = g.get_edge(e_idx);

        std::size_t new_v2{};
        v_type v2_type;

        if (e.get_v2_idx() == v_idx){
          auto [l, r] =  pu::frm_bidirected_idx(e.get_v1_idx());
          new_v2 = e.get_v1_end() ==  v_end::l ? l : r;
          v2_type = e.get_v1_end() == v_end::l ? v_type::l : v_type::r;
        }
        else {
          auto [l, r] =  pu::frm_bidirected_idx(e.get_v2_idx());
          new_v2 = e.get_v2_end() == v_end::l ?  l : r;
          v2_type = e.get_v2_end() == v_end::l ? v_type::l : v_type::r;
        }

        if (added_edges.count(pt::unordered_pair(i_r, new_v2))) { continue; }

        added_edges.insert(pt::unordered_pair(i_r, new_v2));

        if (new_v2 == i_r) {
          std::cout << "i_r "<< new_v2 << std::endl;
        }

        if (new_v2 == i_r) { continue; }

        this->add_edge(i_r, v_type::r, new_v2, v2_type, color::gray);

        //this->get_vertex_mut(i_r).add_edge(this->edges.size() - 1, color::gray);
        //this->get_vertex_mut(new_v2).add_edge(this->edges.size() - 1, color::gray);
      }

    return;
  };

  // add dummy vertex to make it all one SCC
  if (add_dummy_vertices) {
    this->dummy_vertices_.push_back(0);
    this->vertices.push_back(Vertex("d_s", 0, v_type::dummy));
  }

  // add duplicate vertices and black edges
  for (std::size_t i {}; i < g.size(); ++i) {
    auto [idx_l, idx_r] = pu::frm_bidirected_idx(i);
    std::size_t v_id = g.v_idx_to_id(i);
    duplicate_vertices(v_id, idx_l , idx_r);
  }

  // add gray edges between duplicated vertices
  for (std::size_t i {}; i < g.size(); ++i) {
    auto [l, r] =  pu::frm_bidirected_idx(i);
    do_gray_edges(i, g.get_vertex_by_idx(i),  l , r);
  }

  // connect dummy vertices to tips
  if (add_dummy_vertices) {
    for (auto [side, id] : g.tips()) {
      v_type vt = side == v_end::l ? v_type::l : v_type::r;
      auto [l, r] =  pu::frm_bidirected_idx(id);
      std::size_t v1 = vt == v_type::l ? l : r;

      this->add_edge(0, v_type::dummy, v1, vt, color::gray);
    }
  }


  /*
  if (g.tips().size() == 1) {
    // walk down until we find a vertex with more than 2 neighbours
    std::size_t curr_v {};
    std::set<std::size_t> visited;
    visited.insert(0);

    while (this->get_neighbours(curr_v).size() <= 2) {
      for (auto [_, curr_v_] : this->get_neighbours(curr_v)) {
        if (visited.count(curr_v_)) { continue; }
        visited.insert(curr_v_);
        curr_v = curr_v_;
      }
    }

    // more than 2 neighbours, we can't have 0 neighbours
    v_type curr_vt = this->get_vertex(curr_v).get_type();
    this->add_edge(0, v_type::dummy, curr_v, curr_vt, color::gray);
  }
  */

  /*
  if (add_dummy_vertices) {
    this->dummy_vertices_.push_back(this->size());
    this->vertices.push_back(Vertex("d_e", this->size(), v_type::dummy));
e  }


  // connect dummy start to graph starts
  if (add_dummy_vertices) {
    for (auto [side, id] : g.graph_start_nodes()) {
      v_type vt = side == v_end::l ? v_type::l : v_type::r;
      auto [l, r] =  common_fns::frm_bidirected_idx(id);
      std::size_t v1 = vt == v_type::l ? l : r;

     this->add_edge(0, v_type::dummy, v1, vt, color::gray);
    }
  }

  // connect dummy end to graph ends

  if (add_dummy_vertices) {

    std::set<side_n_id_t> orphan_tips = g.get_orphan_tips();
    std::set<side_n_id_t> end_tips = g.graph_end_nodes();

    end_tips.insert(orphan_tips.begin(), orphan_tips.end());


    for (auto [side, id] : end_tips) {
      v_type vt = side == v_end::l ? v_type::l : v_type::r;
      auto [l, r] =  common_fns::frm_bidirected_idx(id);
      std::size_t v1 = vt == v_type::l ? l : r;

            this->add_edge(v1, vt,this->size() - 1, v_type::dummy,  color::gray);

    }
  }

  // connect dummy start to dummy end
  if (add_dummy_vertices) {
    this->add_edge(0, v_type::dummy, this->size() - 1, v_type::dummy, color::gray);
  }
  */

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

std::size_t BVariationGraph::add_edge(std::size_t v1, v_type v1_type,
                                      std::size_t v2, v_type v2_type,
                                      color c) {
  this->edges.push_back(Edge(v1, v1_type, v2, v2_type, c));
  this->get_vertex_mut(v1).add_edge(this->edges.size() - 1, c);
  this->get_vertex_mut(v2).add_edge(this->edges.size() - 1, c);

  return this->edges.size() - 1;
}

std::size_t BVariationGraph::add_edge(std::size_t v1, v_type v1_type,
                                      std::size_t v2, v_type v2_type,
                                      color c, std::string label) {
  this->edges.push_back(Edge(v1, v1_type, v2, v2_type, c, label));
  this->get_vertex_mut(v1).add_edge(this->edges.size() - 1, c);
  this->get_vertex_mut(v2).add_edge(this->edges.size() - 1, c);

  return this->edges.size() - 1;
}

void BVariationGraph::add_vertex(std::string handle_str, v_type vertex_type) {
  vertices.push_back(Vertex(handle_str, vertices.size(), vertex_type));
}

bool BVariationGraph::replace_vertex(std::size_t vertex_idx, std::string handle_str, v_type vertex_type) {
  if (this->get_vertex(vertex_idx).get_handle() != "") {
    return false;
  }

  this->get_vertex_mut(vertex_idx) = Vertex(handle_str, vertex_idx, vertex_type);
  return true;
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
                  this->get_vertex(v_idx).get_handle(),
                  (this->get_vertex(v_idx).get_type() == v_type::l ? "+" : "-") );
  }

  //
  for (std::size_t e_idx{}; e_idx < this->edges.size() ; ++e_idx) {
    auto e = this->edges.at(e_idx);
    if (e.get_color() == color::black) {
      std::cout <<
        std::format("\t{} -- {} [color=\"black\"; label=\"{} {}\"];\n",
                    e.get_v1_idx(), e.get_v2_idx(), e.get_label(),
                    (e.get_eq_class() == pc::UNDEFINED_SIZE_T ? "" : std::to_string(e.get_eq_class())) );
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
                                 t.get_vertex(v).name(),
                                 g.get_vertex(dfs_num_to_vtx[v]).get_handle(),
                                 obe.size(),
                                 ibe.size(),
                                 children.size(), g.get_neighbours(v).size());
      }
    }

    for (auto c : t.get_children(v)) {
      if (t.get_vertex(v).dfs_num() >= t.get_vertex(c).dfs_num()) {
        g.get_vertex( dfs_num_to_vtx[v]).get_handle();
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

    std::size_t v_handle = std::stoull (this->get_vertex(v_idx).get_handle());
    std::size_t n_handle = std::stoull (this->get_vertex(n).get_handle());

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
      t.add_vertex(pst::Vertex{counter, v.get_handle(), v.get_type()});
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
