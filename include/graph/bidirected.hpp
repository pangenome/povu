#ifndef BIDIRECTED_HPP
#define BIDIRECTED_HPP

#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <sys/types.h>
#include <vector>
#include <algorithm>
#include <cstddef>
#include <cstring>
#include <format>
#include <iostream>
#include <iterator>
#include <map>
#include <ostream>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <sys/types.h>
#include <tuple>
#include <unordered_set>
#include <vector>


#include "../common/types.hpp"
#include "../common/utils.hpp"

//TODO: find a proper place for the app config
#include "../../src/cli/app.hpp" // for core::config


namespace povu::bidirected {
namespace pgt = povu::graph_types;
namespace pu = povu::utils;

/**
 *
*/
struct PathInfo {
  std::size_t path_id;
  pgt::or_t strand;
  std::size_t step_index;

  // --------------
  // constructor(s)
  // --------------
  PathInfo(std::size_t path_id, pgt::or_t strand, std::size_t step_index)
    : path_id(path_id),strand(strand), step_index(step_index)  {}
};


// undirected edge
class Edge {
  std::size_t v1_id;
  pgt::v_end v1_end;
  std::size_t v2_id;
  pgt::v_end v2_end;

public:
  // --------------
  // constructor(s)
  // --------------
  Edge(std::size_t v1_id, pgt::v_end v1_end , std::size_t v2_id, pgt::v_end v2_end);

  // ---------
  // getter(s)
  // ---------
  std::size_t get_v1_idx() const;
  pgt::v_end get_v1_end() const;
  std::size_t get_v2_idx() const;
  pgt::v_end get_v2_end() const;
  pgt::side_n_id_t get_other_vertex(std::size_t v_id) const;
};


class Vertex {
  std::string label_; // or sequence
  std::size_t v_id; // this is the sequence name in the GFA file. Maybe Should support strings.

  // indexes to the edge vector in Graph
  std::set<std::size_t> e_l;
  std::set<std::size_t> e_r;

  // references (also colors)
  std::vector<PathInfo> refs_;

public:
  // --------------
  // constructor(s)
  // --------------
  //Vertex(std::size_t v_id);
  Vertex(std::size_t v_id, const std::string& label);

  // ---------
  // getter(s)
  // ---------
  std::size_t id() const;
  const std::string& get_label() const;
  const std::set<std::size_t>& get_edges_l() const;
  const std::set<std::size_t>& get_edges_r() const;

  // ---------
  // setter(s)
  // ---------
  void add_edge_l(std::size_t e_id);
  void add_edge_r(std::size_t e_id);
  void add_ref(std::size_t path_id, pgt::or_t strand, std::size_t step_index);
};


class VariationGraph {
  std::vector<Vertex> vertices;
  std::vector<Edge> edges;
  pu::TwoWayMap<std::size_t, std::size_t> v_id_to_idx_; // TODO: reserve size
  std::map<id_t, std::string> refs_; // a map of ref ids to names
  std::set<pgt::side_n_id_t> tips_;

public:
  // --------------
  // constructor(s)
  // --------------
  //VariationGraph();
  VariationGraph(std::size_t vtx_count, std::size_t edge_count);


  // ---------
  // getter(s)
  // ---------
  std::size_t v_id_to_idx(std::size_t v_id) const;
  std::size_t v_idx_to_id(std::size_t v_idx) const;

  std::size_t size() const; // number of vertices
  std::size_t edge_count() const;
  const std::set<pgt::side_n_id_t>& tips() const;
  const Edge& get_edge(std::size_t e_idx) const;
  const Vertex& get_vertex_by_idx(std::size_t v_idx) const;
  const Vertex& get_vertex_by_id(std::size_t v_id) const;
  Vertex& get_vertex_mut_by_id(std::size_t v_id);


  // ---------
  // setter(s)
  // ---------
  void add_tip(std::size_t v_id, pgt::v_end end);
  void add_vertex(std::size_t v_id, const std::string& label);
  void add_edge(std::size_t v1_id, pgt::v_end v1_end, std::size_t v2_id, pgt::v_end v2_end);
  void add_ref(const std::string &ref_name);

  // other
  void summary() const;
};

typedef VariationGraph VG;

std::vector<povu::bidirected::VG> componetize(const povu::bidirected::VG& g, const core::config& app_config);


} // namespace povu::graph
#endif
