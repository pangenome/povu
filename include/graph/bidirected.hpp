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

#include "../../src/cli/app.hpp" // for core::config TODO: find a proper place for the app config
#include "../common/types.hpp"
#include "../common/utils.hpp"
#include "./spanning_tree.hpp"

namespace povu::bidirected {
namespace pgt = povu::graph_types;
namespace pu = povu::utils;
namespace pst = povu::spanning_tree;
namespace pt = povu::types;
namespace  pc = povu::constants;
using namespace povu::graph_types;
namespace pgt = povu::graph_types;
  //namespace pbd = povu::bidirected;

/**
 *
*/
struct PathInfo {
  pt::id_t path_id;
  pgt::or_t strand;
  pt::idx_t step_index;

  // --------------
  // constructor(s)
  // --------------
  PathInfo(pt::id_t path_id, pgt::or_t strand, pt::id_t step_index)
    : path_id(path_id),strand(strand), step_index(step_index)  {}
};

// undirected edge
// stores the index of the vertex in the graph not the id
class Edge {
  pt::idx_t v1_idx_;
  pgt::v_end v1_end;
  pt::idx_t v2_idx_;
  pgt::v_end v2_end;

public:
  // --------------
  // constructor(s)
  // --------------
  Edge(pt::idx_t v1_idx, pgt::v_end v1_end , pt::idx_t v2_idx, pgt::v_end v2_end);

  // ---------
  // getter(s)
  // ---------
  pt::idx_t get_v1_idx() const;
  pgt::v_end get_v1_end() const;
  pt::idx_t get_v2_idx() const;
  pgt::v_end get_v2_end() const;
  pgt::side_n_id_t get_other_vertex(pt::idx_t v_id) const;
};


class Vertex {
  pt::id_t v_id; // this is the sequence name in the GFA file. Maybe Should support strings.
  std::string label_; // or sequence

  // indexes to the edge vector in Graph
  std::set<pt::idx_t> e_l;
  std::set<pt::idx_t> e_r;

  // references (also colors)
  std::vector<PathInfo> refs_;

public:
  // --------------
  // constructor(s)
  // --------------
  Vertex(pt::id_t v_id, const std::string& label = "");

  // ---------
  // getter(s)
  // ---------
  pt::id_t id() const;
  const std::string& get_label() const;
  const std::set<pt::idx_t>& get_edges_l() const;
  const std::set<pt::idx_t>& get_edges_r() const;

  // ---------
  // setter(s)
  // ---------
  void add_edge_l(pt::idx_t e_idx);
  void add_edge_r(pt::idx_t e_idx);
  void add_ref(pt::id_t path_id, pgt::or_t strand, pt::idx_t step_index);
};


class VariationGraph {
  std::vector<Vertex> vertices;
  std::vector<Edge> edges;
  pu::TwoWayMap<std::size_t, std::size_t> v_id_to_idx_; // TODO: reserve size
  std::map<id_t, std::string> refs_; // a map of ref ids to names
  std::set<pgt::side_n_id_t> tips_; // the set of side and id of the tips

  pt::idx_t dummy_idx_; // index of the dummy vertex

public:
  // --------------
  // constructor(s)
  // --------------
  VariationGraph(pt::idx_t vtx_count, pt::idx_t edge_count);


  // ---------
  // getter(s)
  // ---------
  pt::idx_t v_id_to_idx(pt::id_t v_id) const;
  pt::id_t v_idx_to_id(pt::idx_t v_idx) const;

  pt::idx_t vtx_count() const;
  pt::idx_t edge_count() const;
  const std::set<pgt::side_n_id_t>& tips() const;
  const Edge& get_edge(pt::idx_t e_idx) const;
  const Vertex& get_vertex_by_idx(pt::idx_t v_idx) const;
  const Vertex& get_vertex_by_id(pt::id_t v_id) const;
  Vertex& get_vertex_mut_by_id(pt::id_t v_id);

  // TODO: reduce number of dummy related methods
  const Vertex& get_dummy_vertex() const;
  std::size_t get_dummy_idx() const;
  bool has_dummy() const;


  // ---------
  // setter(s)
  // ---------
  void add_tip(std::size_t v_id, pgt::v_end end); // TODO: give a name that shows that we don't add a vertex but mark a vertex as tip
  // returns the index (v_idx) of the added vertex
  pt::idx_t add_vertex(pt::id_t v_id, const std::string& label);

  pt::idx_t add_edge(pt::id_t v1_id, pgt::v_end v1_end, pt::id_t v2_id, pgt::v_end v2_end);
  void add_ref(const std::string &ref_name);
  /**
   * make sure not to run more than once
   *
   *
   */
  void untip(); //
  void shrink_to_fit();

  // other
  void summary() const;
  void print_dot(std::ostream& os) const;
};

typedef VariationGraph VG;

std::vector<VG *> componetize(const VG &g);

pst::Tree compute_spanning_tree(const VG &g);

} // namespace povu::bidirected
#endif
