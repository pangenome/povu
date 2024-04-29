#include <cstddef>
#include <format>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <vector>
#include <stack>

#include "./biedged.hpp"


namespace biedged {

using namespace graph_types;

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
  black_edge(constants::UNDEFINED_SIZE_T),
  grey_edges(std::set<std::size_t>()),
  type(v_type::l),
  handle(std::string()),
  vertex_idx(0)
{}

Vertex::Vertex(const std::string& id, std::size_t vertex_idx, v_type vertex_type) :
  black_edge(constants::UNDEFINED_SIZE_T),
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

std::set<std::size_t> Vertex::get_neighbours() const {
  std::set<std::size_t> neighbours { this->get_grey_edges().begin(), this->get_grey_edges().end() };
  if (this->get_black_edge() != constants::UNDEFINED_SIZE_T) {
    neighbours.insert(this->get_black_edge());
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
      if (bg.get_vertex(i).get_black_edge() == constants::UNDEFINED_SIZE_T) {
        std::cout << "no black edge failed check i " << i << std::endl;
      }
    }
  }
}

BVariationGraph::BVariationGraph(const bidirected::VariationGraph &g, bool add_dummy_vertices) {
  std::string fn_name = std::format("[povu::BVariationGraph::{}]", __func__);

  if (add_dummy_vertices) {
    if (g.graph_start_nodes().empty()) {
      throw std::domain_error(fn_name + " : graph has no start nodes");
    }

    if (g.graph_end_nodes().empty()) {
      throw std::domain_error(fn_name + " : graph has no end nodes");
    }
  }

  std::size_t biedged_size = g.size() * 2 + (add_dummy_vertices ? 2 : 0);

  this->vertices = std::vector<biedged::Vertex>();
  this->vertices.reserve(biedged_size);

  // this->vertices = std::vector<biedged::Vertex>(biedged_size, biedged::Vertex());

  auto new_l = [&](std::size_t i) { return (2*(i+1)) - 2 + (add_dummy_vertices ?  1 : 0 ); };
  auto new_r = [&](std::size_t i) { return (2*(i+1)) - 1 + (add_dummy_vertices ?  1 : 0 ); };

  auto duplicate_vertices = [this, g, add_dummy_vertices](bidirected::Vertex const& current_vg_vertex, std::size_t i_l, std::size_t i_r) {

    // --------------------
    // add vertices + and -
    // --------------------

    {
      this->vertices.push_back(Vertex(current_vg_vertex.get_name(), i_l, v_type::l));
      this->vertices.push_back(Vertex(current_vg_vertex.get_name(), i_r, v_type::r));
    }

    // ----------------
    // add a black edge
    // ----------------
    {
      std::size_t e_idx =
        this->add_edge(i_l, v_type::l,
                       i_r, v_type::r,
                       color::black,
                       current_vg_vertex.get_label());

      this->get_vertex_mut(i_l).add_edge(e_idx, color::black);
      this->get_vertex_mut(i_r).add_edge(e_idx, color::black);
    }
  };

  std::set<utils::unordered_pair> added_edges;

  /*
    add gray edges
  */
  auto do_gray_edges = [&](std::size_t v_idx, bidirected::Vertex const& current_vg_vertex, std::size_t i_l, std::size_t i_r) {
      // add gray edges incident with the 5' (left/+) vertex
      for (auto e_idx : current_vg_vertex.get_edges_l()) {
        const bidirected::Edge &e = g.get_edge(e_idx);

        std::size_t new_v1{};
        v_type v1_type{};

        if (e.get_v1_idx() == v_idx) {
          auto [l, r] =  common_fns::frm_bidirected_idx(e.get_v2_idx());
          new_v1 = e.get_v2_end() == v_end::r ? r : l;
          v1_type = e.get_v2_end() == v_end::l ? v_type::r : v_type::l;
        }
        else {
          auto [l, r] =  common_fns::frm_bidirected_idx(e.get_v1_idx());
          new_v1 = e.get_v1_end() == v_end::l ? l : r;
          v1_type = e.get_v1_end() == v_end::l ? v_type::l : v_type::r;
        }

        if (added_edges.count(utils::unordered_pair(new_v1, i_l))) { continue; }

        added_edges.insert(utils::unordered_pair(new_v1, i_l));

        if (new_v1 == i_l) {
          std::cout << "e " << e << std::endl;
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

        bidirected::Edge e = g.get_edge(e_idx);

        std::size_t new_v2{};
        v_type v2_type;

        if (e.get_v2_idx() == v_idx){
          auto [l, r] =  common_fns::frm_bidirected_idx(e.get_v1_idx());
          new_v2 = e.get_v1_end() ==  v_end::l ? l : r;
          v2_type = e.get_v1_end() == v_end::l ? v_type::l : v_type::r;
        }
        else {
          auto [l, r] =  common_fns::frm_bidirected_idx(e.get_v2_idx());
          new_v2 = e.get_v2_end() == v_end::l ?  l : r;
          v2_type = e.get_v2_end() == v_end::l ? v_type::l : v_type::r;
        }

        if (added_edges.count(utils::unordered_pair(i_r, new_v2))) { continue; }

        added_edges.insert(utils::unordered_pair(i_r, new_v2));

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
  for (std::size_t i = 0; i < g.size(); ++i) {
    auto [l, r] =  common_fns::frm_bidirected_idx(i);
    duplicate_vertices(g.get_vertex(i), l , r);
  }

  // add gray edges between duplicated vertices
  for (std::size_t i = 0; i < g.size(); ++i) {
    auto [l, r] =  common_fns::frm_bidirected_idx(i);
    do_gray_edges(i, g.get_vertex(i),  l , r);
  }

  // connect dummy vertices to tips
  if (add_dummy_vertices) {
    for (auto [side, id] : g.tips()) {
      v_type vt = side == v_end::l ? v_type::l : v_type::r;
      auto [l, r] =  common_fns::frm_bidirected_idx(id);
      std::size_t v1 = vt == v_type::l ? l : r;

      this->add_edge(0, v_type::dummy, v1, vt, color::gray);
    }
  }

  /*
  if (add_dummy_vertices) {
    this->dummy_vertices_.push_back(this->size());
    this->vertices.push_back(Vertex("d_e", this->size(), v_type::dummy));
  }


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

std::set<std::pair<color, std::size_t>> BVariationGraph::get_neighbours(std::size_t vertex_idx) const {
  std::set<std::pair<color, std::size_t>> neighbours;
  Vertex const& v = this->get_vertex(vertex_idx);

  for (auto e_idx : v.get_grey_edges()) {
    const Edge& e = this->get_edge(e_idx);
    std::size_t other_vertex = e.get_other_vertex(vertex_idx);
    std::pair<color, std::size_t> p = std::make_pair(color::gray, other_vertex);
    neighbours.insert(p);
  }

  std::size_t black_e_idx = v.get_black_edge();
  if (black_e_idx != constants::UNDEFINED_SIZE_T) {
  neighbours.insert(std::make_pair(
      color::black,
      this->get_edge(black_e_idx).get_other_vertex(vertex_idx)));
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

void BVariationGraph::update_eq_classes(spanning_tree::Tree& st) {
  const std::map<std::size_t, std::pair<spanning_tree::EdgeType, std::size_t>>& m = st.get_g_edge_idx_map();

  for (auto& [e_idx, v] : m) {
    auto [e_type, t_e_idx] = v;

    std::size_t eq_class = (e_type == spanning_tree::EdgeType::tree_edge)
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
                    (e.get_eq_class() == constants::UNDEFINED_SIZE_T ? "" : std::to_string(e.get_eq_class())) );
    }
    else {
      std::cout <<
        std::format("\t{} -- {} [color=\"gray\" label=\"{}\"];\n", e.get_v1_idx(), e.get_v2_idx(),
                    (e.get_eq_class() == constants::UNDEFINED_SIZE_T ? "" : std::to_string(e.get_eq_class())));
    }
  }

  std::cout << "}" << std::endl;
}

// TODO: make ref params const
void validate_spanning_tree(BVariationGraph const& g,
                            spanning_tree::Tree& t,
                            std::map<std::size_t, std::size_t>& dfs_num_to_vtx) {
  if (t.size() != g.size()) {
    std::cerr << "size mismatch " << t.size() << " " << g.size() << "\n";
  }

  // compare edge count
  if (t.tree_edge_count() + t.back_edge_count() != g.num_edges()) {
    std::cerr << "edge count mismatch " << t.tree_edge_count() << " " << t.back_edge_count() << " " << g.num_edges() << "\n";
  }

  for (std::size_t v{}; v < t.size() ; v++) {
    std::set<std::size_t>const& children = t.get_children(v);
    std::set<std::size_t>const& obe = t.get_obe(v);
    std::set<std::size_t>const& ibe = t.get_ibe(v);

    if (t.is_root(v)) {
      if (obe.size() + ibe.size() + children.size() != g.get_neighbours(dfs_num_to_vtx[v]).size()) {
        std::cerr << "neighbour mismatch (root) " << v << " " << obe.size() << " " << ibe.size() << " " << children.size() << " " << g.get_neighbours(v).size() << "\n";
      }
    }
    else {
      if (obe.size() + ibe.size() + children.size() + 1 != g.get_neighbours(dfs_num_to_vtx[v]).size()) {
        std::cerr << "neighbour mismatch " << v << " " << obe.size() << " " << ibe.size() << " " << children.size() << " " << g.get_neighbours(v).size() << "\n";
      }
    }

    for (auto c : t.get_children(v)) {
      if (t.get_vertex(v).dfs_num() >= t.get_vertex(c).dfs_num()) {
        std::cerr << "te weird " << v << " " << t.get_vertex(v).dfs_num() << " "
                  << c << " " << t.get_vertex(c).dfs_num() << "\n";
      }
    }

    for (std::size_t tgt : t.get_obe(v)) {
      if (t.get_vertex(v).dfs_num() <= t.get_vertex(tgt).dfs_num()) {
        std::cerr << "obe weird " << v << " " << t.get_vertex(v).dfs_num() << " "
                  << tgt << " " << t.get_vertex(tgt).dfs_num() << "\n";
      }
    }

    for (std::size_t src : t.get_ibe(v)) {
      if (t.get_vertex(v).dfs_num() >= t.get_vertex(src).dfs_num()) {
        std::cerr << "ibe weird" << v << " " << t.get_vertex(v).dfs_num() << " "
                  << src << " " << t.get_vertex(src).dfs_num() << "\n";
      }
    }
  }
}


spanning_tree::Tree BVariationGraph::compute_spanning_tree() const {
  spanning_tree::Tree t = spanning_tree::Tree(this->size());

  std::stack<std::tuple<std::size_t, std::size_t, bidirected::color>> s; // parent and v idx
  std::set<std::size_t> visited;
  std::size_t v_idx {}; // set start node to 0
  s.push({constants::INVALID_ID, v_idx, bidirected::color::gray});
  visited.insert(v_idx);
  std::size_t counter{}; // dfs pre visit counter

  std::map<std::size_t, std::size_t> vtx_to_dfs_num;
  std::map<std::size_t, std::size_t> dfs_num_to_vtx;
  std::set<std::size_t> in_tree; // graph vertices in the tree
  std::set<utils::unordered_pair> back_edges; // avoids duplicate back edges

  auto is_parent = [&](std::size_t p, std::size_t c) -> bool {
    return t.get_children(vtx_to_dfs_num[p]).count(vtx_to_dfs_num[c]);
  };

  auto is_parent_child = [&](std::size_t n, std::size_t v_idx) -> bool {
    return is_parent(n, v_idx) || is_parent(v_idx, n);
  };

  while (!s.empty()) {
    auto [p_idx, v_idx, c] = s.top();

    biedged::Vertex const& v = this->get_vertex(v_idx);

    if (in_tree.find(v_idx) == in_tree.end()) {
      t.add_vertex(spanning_tree::Vertex{counter, v.get_handle(), v.get_type()});
      in_tree.insert(v_idx);
      vtx_to_dfs_num[v_idx] = counter;
      dfs_num_to_vtx[counter] = v_idx;

      if (p_idx != constants::INVALID_ID) {
        t.add_tree_edge(vtx_to_dfs_num[p_idx], counter, constants::UNDEFINED_SIZE_T, c);
      }

      counter++;
    }

    bool explored {true};
    for (auto [c, n] : this->get_neighbours(v_idx)) {
      if (visited.find(n) == visited.end()) {
        s.push({v_idx, n, c});
        visited.insert(n);
        explored = false;
        break;
      }
      else if (!is_parent_child(n, v_idx) && back_edges.find({n, v_idx}) == back_edges.end())  {
        t.add_be(vtx_to_dfs_num[v_idx], vtx_to_dfs_num[n],constants::UNDEFINED_SIZE_T, false, c);
        back_edges.insert({v_idx, n});
      }
    }

    if (explored) { s.pop(); }
  }

#ifdef DEBUG   // check the correctness of the tree
  validate_spanning_tree(*this, t, dfs_num_to_vtx);
#endif
  return t;
}

}; // namespace biedged
