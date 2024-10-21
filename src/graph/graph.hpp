#ifndef PV_GRAPH_HPP
#define PV_GRAPH_HPP
#include <tuple>
#include <vector>
#include <set>
#include <map>

#include "../common/types.hpp"
#include "../common/utils.hpp"
#include "../cli/app.hpp"

namespace povu::graph {
namespace pgt = povu::graph_types;
namespace pu = povu::utils;

// undirected edge
class Edge {
  std::size_t v1_id;
  pgt::v_end v1_end;
  std::size_t v2_id;
  pgt::v_end v2_end;

public:
  Edge(std::size_t v1_id, pgt::v_end v1_end , std::size_t v2_id, pgt::v_end v2_end);
  std::size_t get_v1_idx() const;
  pgt::v_end get_v1_end() const;
  std::size_t get_v2_idx() const;
  pgt::v_end get_v2_end() const;

  pgt::side_n_id_t get_other_vertex(std::size_t v_id) const;
};


class Vertex {
  std::size_t v_id;
  std::set<std::size_t> e_l;
  std::set<std::size_t> e_r;

public:
  Vertex(std::size_t v_id);
  std::size_t id() const;
  void add_edge_l(std::size_t e_id);
  void add_edge_r(std::size_t e_id);
  const std::set<std::size_t>& get_edges_l() const;
  const std::set<std::size_t>& get_edges_r() const;
};

class Graph {
  std::vector<Vertex> vertices;
  std::vector<Edge> edges;
  pu::TwoWayMap<std::size_t, std::size_t> v_id_to_idx_; // TODO: reserve size
  std::set<pgt::side_n_id_t> tips_;

public:
  // constructors
  Graph();
  Graph(std::size_t v_count, std::size_t e_count);

  // getters

  std::size_t v_id_to_idx(std::size_t v_id) const;
  std::size_t v_idx_to_id(std::size_t v_idx) const;

  std::size_t size() const; // number of vertices
  std::size_t edge_count() const;
  const std::set<pgt::side_n_id_t>& tips() const;
  const Edge& get_edge(std::size_t e_idx) const;
  const Vertex& get_vertex_by_idx(std::size_t v_id) const;
  const Vertex& get_vertex_by_id(std::size_t v_idx) const;

  // setters
  void add_tip(std::size_t v_id, pgt::v_end end);
  void add_vertex(std::size_t v_id);
  void add_edge(std::size_t v1_id, pgt::v_end v1_end, std::size_t v2_id, pgt::v_end v2_end);

  // other
  void summary() const;
};

std::vector<povu::graph::Graph> componetize(const povu::graph::Graph& g, const core::config& app_config);

} // namespace povu::graph
#endif
