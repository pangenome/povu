#include <cstddef>
#include <format>
#include <iostream>
#include <iterator>
#include <string>
#include <vector>
#include <stack>

#include "./biedged.hpp"
#include "./bidirected.hpp"

namespace biedged {

/*
  VertexType
 */

// implement << operator for VertexType
std::ostream& operator<<(std::ostream& os, const VertexType& vt) {
    switch (vt) {
    case VertexType::l:
    os << "+";
    break;
    case VertexType::r:
    os << "-";
    break;
    default:
    os << "*";
    break;
    }

    return os;
}

/*
  Edge
*/
Edge::Edge(std::size_t v1, VertexType v1_type,
           std::size_t v2, VertexType v2_type,
           core::color c)
  : v1_idx(v1), v1_type(v1_type), v2_idx(v2), v2_type(v2_type), c(c)
{
  std::string fn_name = "[povu::graph::biedged::Edge]";

  if (v1 == v2) {
    throw std::invalid_argument(std::format("{} Self-loops are not allowed {} {}", fn_name, v1, v2));
  }

  if (c == core::color::black) {
    // throw an argument error if the edge is black and no label is provided
    throw std::invalid_argument(fn_name + " edges must have a label");
  }

  this->label = std::string();
}

Edge::Edge(std::size_t v1, VertexType v1_type,
           std::size_t v2, VertexType v2_type,
           core::color c, std::string label)
  : v1_idx(v1), v1_type(v1_type), v2_idx(v2), v2_type(v2_type), c(c), label(label)
{
  std::string fn_name = "[povu::graph::biedged::Edge]";

  if (v1 == v2) {
    throw std::invalid_argument(fn_name + " Self-loops are not allowed");
  }

  if (c == core::color::gray) {
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

core::color Edge::get_color() const {
  return this->c;
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
    throw std::invalid_argument(
                                std::format("{} Vertex {} is not part of edge", fn_name, vertex_index));
  }
}

const std::string& Edge::get_label() const {
    return this->label;
}

void Edge::set_v1_idx(std::size_t i) {
    this->v1_idx = i;
}

void Edge::set_v2_idx(std::size_t i) {
    this->v2_idx = i;
}

std::ostream& operator<<(std::ostream& os, const Edge& e) {
    os << "Edge(" << e.v1_idx << ", "
       << (e.v1_type == VertexType::l ? "-" : "+")
       << e.v2_idx << ", "
       << (e.v2_type == VertexType::l ? "-" : "+")
       << e.c << ")";
    return os;
}

/*
  Vertex
*/
Vertex::Vertex() :
  black_edge(core::constants::UNDEFINED_SIZE_T),
  grey_edges(std::set<std::size_t>()),
  paths(std::vector<PathInfo>()),
  type(VertexType::l),
  handle(std::string()),
  is_reversed_(false),
  vertex_idx(0)
{}

Vertex::Vertex(const std::string& id, std::size_t vertex_idx, VertexType vertex_type) :
  black_edge(core::constants::UNDEFINED_SIZE_T),
  grey_edges(std::set<std::size_t>()),
  paths(std::vector<PathInfo>()),
  type(vertex_type),
  handle(id),
  is_reversed_(false),
  vertex_idx(vertex_idx)
{}


const std::string& Vertex::get_handle() const {
    return this->handle;
}

std::set<std::size_t> Vertex::get_grey_edges() const {
    return this->grey_edges;
}

std::size_t Vertex::get_black_edge() const {
    return this->black_edge;
}

std::size_t Vertex::get_vertex_idx() const {
  return this->vertex_idx;
}

void Vertex::add_edge(std::size_t edge_idx, core::color c) {
  if (c == core::color::black) {
    black_edge = edge_idx;
  } else {
    grey_edges.insert(edge_idx);
  }
}

VertexType Vertex::get_type() const {
    return type;
}

void Vertex::set_vertex_idx(std::size_t i) {
    this->vertex_idx = i;
}

/*
  BVariationGraph
 */

std::size_t BVariationGraph::add_edge(std::size_t v1, VertexType v1_type,
                                      std::size_t v2, VertexType v2_type,
                                      core::color c) {
  this->edges.push_back(Edge(v1, v1_type, v2, v2_type, c));
  this->get_vertex_mut(v1).add_edge(this->edges.size() - 1, c);
  this->get_vertex_mut(v2).add_edge(this->edges.size() - 1, c);

  return this->edges.size() - 1;
}

std::size_t BVariationGraph::add_edge(std::size_t v1, VertexType v1_type,
                                      std::size_t v2, VertexType v2_type,
                                      core::color c, std::string label) {
  this->edges.push_back(Edge(v1, v1_type, v2, v2_type, c, label));
  this->get_vertex_mut(v1).add_edge(this->edges.size() - 1, c);
  this->get_vertex_mut(v2).add_edge(this->edges.size() - 1, c);

  return this->edges.size() - 1;
}

void BVariationGraph::add_vertex(std::string handle_str, VertexType vertex_type) {
  vertices.push_back(Vertex(handle_str, vertices.size(), vertex_type));
}

bool BVariationGraph::replace_vertex(std::size_t vertex_idx, std::string handle_str, VertexType vertex_type) {
    if (this->get_vertex(vertex_idx).get_handle() != "") {
        return false;
    }

    this->get_vertex_mut(vertex_idx) = Vertex(handle_str, vertex_idx, vertex_type);
    return true;
}

const biedged::Vertex& BVariationGraph::get_vertex(std::size_t i) const {
    return vertices.at(i);
}

std::size_t BVariationGraph::size() const {
    return vertices.size();
}

Vertex& BVariationGraph::get_vertex_mut(std::size_t i) {
    return vertices.at(i);
}

BVariationGraph::BVariationGraph(const bidirected::VariationGraph &g, bool add_dummy_vertices) {

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
      this->vertices.push_back(Vertex(current_vg_vertex.get_name(), i_l, VertexType::l));
      this->vertices.push_back(Vertex(current_vg_vertex.get_name(), i_r, VertexType::r));
    }

    // ----------------
    // add a black edge
    // ----------------
    {
      std::size_t e_idx =
        this->add_edge(i_l, VertexType::l,
                       i_r, VertexType::r,
                       core::color::black,
                       current_vg_vertex.get_label());

      this->get_vertex_mut(i_l).add_edge(e_idx, core::color::black);
      this->get_vertex_mut(i_r).add_edge(e_idx, core::color::black);
    }
  };

  std::set<unordered_pair> added_edges;

  /*
    add gray edges
  */
  auto do_gray_edges = [&](std::size_t v_idx, bidirected::Vertex const& current_vg_vertex, std::size_t i_l, std::size_t i_r) {
      // add gray edges incident with the 5' (left/+) vertex
      for (auto e_idx : current_vg_vertex.get_edges_l()) {
        const bidirected::Edge &e = g.get_edge(e_idx);

        std::size_t new_v1{};
        VertexType v1_type{};

        if (e.get_v1_idx() == v_idx) {
          new_v1 = e.get_v2_end() == bidirected::VertexEnd::r ? new_r(e.get_v2_idx()) : new_l(e.get_v2_idx());
          // std::cout << "neq " << new_v1 << "\n";
          v1_type = e.get_v2_end() == bidirected::VertexEnd::l ? VertexType::r : VertexType::l;
        } else {
          new_v1 = e.get_v1_end() == bidirected::VertexEnd::l ? new_l(e.get_v1_idx()) : new_r(e.get_v1_idx());
          // std::cout << "nneq " << new_v1 << "\n";
          v1_type = e.get_v1_end() == bidirected::VertexEnd::l ? VertexType::l : VertexType::r;
        }

        if (added_edges.count(unordered_pair(new_v1, i_l))) { continue; }

        added_edges.insert(unordered_pair(new_v1, i_l));

        if (new_v1 == i_l) {
          std::cout << "e " << e << std::endl;
          std::cout << "v_idx " << v_idx << " " << e.get_v1_idx() << " " << e.get_v2_end() << std::endl;
          std::cout << "n " << new_v1 << "r " << new_r(e.get_v2_idx()) << " l " << new_l(e.get_v2_idx()) << std::endl;
          //continue;
        }

        std::size_t gray_e_idx = this->add_edge(new_v1, v1_type, i_l, VertexType::l,core::color::gray);

        //this->get_vertex_mut(new_v1).add_edge(gray_e_idx, core::color::gray);
        //this->get_vertex_mut(i_l).add_edge(gray_e_idx, core::color::gray);
      }

      // add gray edges incident with the 3' (right/-) vertex
      for (auto e_idx : current_vg_vertex.get_edges_r()) {

        bidirected::Edge e = g.get_edge(e_idx);

        std::size_t new_v2{};
        VertexType v2_type;

        if (e.get_v2_idx() == v_idx){
          new_v2 = e.get_v1_end() == bidirected::VertexEnd::l ? new_l(e.get_v1_idx()) : new_r(e.get_v1_idx());
          v2_type = e.get_v1_end() == bidirected::VertexEnd::l ? VertexType::l : VertexType::r;
        }
        else {
          new_v2 =
            e.get_v2_end() == bidirected::VertexEnd::l ?
            new_l(e.get_v2_idx()) :
            new_r(e.get_v2_idx());

          v2_type =
            e.get_v2_end() == bidirected::VertexEnd::l ? VertexType::l : VertexType::r;
        }

        if (added_edges.count(unordered_pair(i_r, new_v2))) { continue; }

        added_edges.insert(unordered_pair(i_r, new_v2));

                if (new_v2 == i_r) {
                  std::cout << "i_r "<< new_v2 << std::endl;
        }

                if (new_v2 == i_r) {
                  continue;
                }

        this->add_edge(i_r, VertexType::r,
                       new_v2, v2_type,
                       core::color::gray);

        //this->get_vertex_mut(i_r).add_edge(this->edges.size() - 1, core::color::gray);
        //this->get_vertex_mut(new_v2).add_edge(this->edges.size() - 1, core::color::gray);
      }

    return;
  };

  if (add_dummy_vertices) { // add dummy start
    this->vertices.push_back(Vertex("d_s", 0, VertexType::dummy));
  }

  for (std::size_t i = 0; i < g.size(); ++i) {
    duplicate_vertices(g.get_vertex(i), new_l(i) , new_r(i));
  }

  if (add_dummy_vertices) { // add grey edges from dummy start to graph starts
    for (auto i : g.graph_start_nodes()) {
      VertexType vt = g.get_vertex(i).get_edges_l().empty() ? VertexType::l : VertexType::r ;
      std::size_t v2 = vt == VertexType::l ? new_l(i) : new_r(i);

      this->add_edge(0, VertexType::dummy, v2, vt, core::color::gray);
    }
  }

  for (std::size_t i = 0; i < g.size(); ++i) {
    do_gray_edges(i, g.get_vertex(i),  new_l(i) , new_r(i));
  }

  if (add_dummy_vertices) {
    // add dummy end
    this->vertices.push_back(Vertex("d_e", this->size(), VertexType::dummy));

    // add grey edges from graph ends to dummy end
    for (auto i : g.graph_end_nodes()) {
      VertexType vt = g.get_vertex(i).get_edges_r().empty() ? VertexType::r : VertexType::l ;
      std::size_t v1 = vt == VertexType::l ? new_l(i) : new_r(i);
      this->add_edge(v1, vt, this->size() - 1, VertexType::dummy, core::color::gray);
    }

    // add an egde between the dummy start and the dummy end
    if (add_dummy_vertices) {
      this->add_edge(0, VertexType::dummy, this->size() - 1, VertexType::dummy, core::color::gray);
    }
  }
}

void BVariationGraph::print_dot() const {
  std::cout << "graph G {\n" <<
    "\trankdir=LR;\n" <<
    "\tnode [shape=circle];\n";

  for (std::size_t v_idx{}; v_idx < this->vertices.size(); ++v_idx) {
    std::cout <<
      std::format("\t{} [label=\"{}{}\"]\n",
                  v_idx,
                  this->get_vertex(v_idx).get_handle(),
                  (this->get_vertex(v_idx).get_type() == VertexType::l ? "+" : "-") );
  }

  //
  for (std::size_t e_idx{}; e_idx < this->edges.size() ; ++e_idx) {
    auto e = this->edges.at(e_idx);
    if (e.get_color() == core::color::black) {
      std::cout <<
        std::format("\t{} -- {} [color=\"black\"; label=\"{}\"];\n",
                    e.get_v1_idx(), e.get_v2_idx(), e.get_label());
    }
    else {
      std::cout <<
        std::format("\t{} -- {} [color=\"gray\"];\n", e.get_v1_idx(), e.get_v2_idx());
    }
  }

  std::cout << "}" << std::endl;
}

spanning_tree::Tree BVariationGraph::compute_spanning_tree() const {
  std::string fn_name = std::format("[povu::biedged::{}]", __func__);

  spanning_tree::Tree t = spanning_tree::Tree(this->size());

  std::size_t v_idx {}; // set start node to 0
  std::stack<std::size_t> s;
  s.push(v_idx);

  std::set<std::size_t> visited;
  std::size_t counter {};

  while (!s.empty()) {
    v_idx = s.top();

    if (visited.find(v_idx) == visited.end()) {
      t.set_dfs_num(v_idx, counter);
      t.set_sort(counter, v_idx);
      t.set_sort_g(v_idx, counter);
      t.get_vertex_mut(v_idx).set_name(this->get_vertex(v_idx).get_handle()) ;
      ++counter;
    }

    visited.insert(v_idx);

    bool explored{true};
    Vertex const &v = this->get_vertex(v_idx);

    std::set<size_t> adj_edges = v.get_grey_edges();
    // if (v.get_black_edge() != core::constants::UNDEFINED_SIZE_T) { adj_edges.insert(v.get_black_edge()); }
    if (v.get_type() != VertexType::dummy) { adj_edges.insert(v.get_black_edge()); }

    for (auto e_idx : adj_edges) {
      const Edge &e = this->edges.at(e_idx);
      std::size_t adj_v_idx = e.get_other_vertex(v_idx);

      if (visited.find(adj_v_idx) == visited.end()) {
        t.add_tree_edge(v_idx, adj_v_idx, e.get_color());
        s.push(adj_v_idx);
        explored = false;
        break;
      }
      else if (
        !t.is_root(v_idx) &&
        t.get_parent(v_idx) != adj_v_idx &&
        !t.has_child(v_idx, adj_v_idx) &&
        !t.has_ibe(v_idx, adj_v_idx) &&
        !t.has_obe(v_idx, adj_v_idx)
      ) {
        // TODO: why the has child and not parent test?
        //std::cout << "adding back edge: " << v_idx << " -> " << a << std::endl;
        t.add_be(v_idx, adj_v_idx, false, e.get_color());
      }
    }

    if (explored) { s.pop(); }
  }

  for (std::size_t i{}; i < this->size(); ++i) {
    if (t.get_vertex(i).is_null()) {
      throw std::logic_error(std::format("{}: vertex {} is null. {}/{} vertices explored\n",
                                         fn_name, i, counter, this->size()));
    }
  }

  return t;
}

}; // namespace biedged
