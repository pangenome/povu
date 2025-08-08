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
#include "../common/utils.hpp"
#include "./spanning_tree.hpp"

namespace povu::bidirected {
inline constexpr std::string_view MODULE = "povu::bidirected";

namespace pu = povu::utils;
namespace pst = povu::spanning_tree;
namespace pt = povu::types;
namespace pc = povu::constants;
using namespace povu::types::graph;
namespace pgt = povu::types::graph;


/**
 *
*/
struct PathInfo {
  pt::id_t path_id;
  pgt::or_e strand;
  pt::idx_t step_index;

  // --------------
  // constructor(s)
  // --------------
  PathInfo(pt::id_t path_id, pgt::or_e strand, pt::id_t step_index)
    : path_id(path_id),strand(strand), step_index(step_index)  {}
};
typedef PathInfo pi_t;

typedef PathInfo VtxRefInfo; 

// undirected edge
// stores the index of the vertex in the graph not the id
class Edge {
  pt::idx_t v1_idx_;
  pgt::v_end_e v1_end;
  pt::idx_t v2_idx_;
  pgt::v_end_e v2_end;

public:
  // --------------
  // constructor(s)
  // --------------
  Edge(pt::idx_t v1_idx, pgt::v_end_e v1_end , pt::idx_t v2_idx, pgt::v_end_e v2_end);

  // ---------
  // getter(s)
  // ---------
  pt::idx_t get_v1_idx() const;
  pt::idx_t &get_v1_idx_mut();
  pgt::v_end_e get_v1_end() const;
  pt::idx_t get_v2_idx() const;
  pt::idx_t &get_v2_idx_mut();
  pgt::v_end_e get_v2_end() const;
  pgt::side_n_idx_t get_other_vtx(pt::idx_t v_id, pgt::v_end_e v_end) const;
  pgt::side_n_idx_t get_other_vtx(pt::idx_t v_id) const; // if you don't care for self loops
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
  std::string get_rc_label() const; // reverse complement of the label
  const std::set<pt::idx_t>& get_edges_l() const;
  const std::set<pt::idx_t>& get_edges_r() const;
  const std::vector<PathInfo>& get_refs() const;

  // ---------
  // setter(s)
  // ---------
  void add_edge_l(pt::idx_t e_idx);
  void add_edge_r(pt::idx_t e_idx);
  void add_ref(pt::id_t ref_id, pgt::or_e strand, pt::idx_t step_index);
};


class VariationGraph {
  std::vector<Vertex> vertices;
  std::vector<Edge> edges;
  pu::TwoWayMap<std::size_t, std::size_t> v_id_to_idx_; // TODO: reserve size
  //std::map<id_t, std::string> refs_; // a map of ref ids to names
  //std::map<id_t, pgt::RefInfo> ref_info_; // a map of ref ids to ref info
  bool has_refs_;
  pgt::Refs refs_;
  std::set<pgt::side_n_id_t> tips_; // the set of side and id of the tips

public:
  // --------------
  // constructor(s)
  // --------------
  VariationGraph(pt::idx_t vtx_count, pt::idx_t edge_count, bool inc_refs);


  // ---------
  // getter(s)
  // ---------
  pt::idx_t v_id_to_idx(pt::id_t v_id) const;
  pt::id_t v_idx_to_id(pt::idx_t v_idx) const;

  pt::idx_t vtx_count() const;
  pt::idx_t edge_count() const;
  const std::set<pgt::side_n_id_t>& tips() const;
  const Edge& get_edge(pt::idx_t e_idx) const;
  Edge & get_edge_mut(pt::idx_t e_idx);
  // TODO replace vertex with v?
  const Vertex& get_vertex_by_idx(pt::idx_t v_idx) const;
  const Vertex& get_vertex_by_id(pt::id_t v_id) const;
  Vertex& get_vertex_mut_by_id(pt::id_t v_id);

  // ref
  const std::string &get_ref_label(pt::id_t ref_id) const;
  const pgt::Ref &get_ref_by_id(pt::id_t ref_id) const;
  pgt::Ref &get_ref_by_id_mut(pt::id_t ref_id);
  pt::id_t get_ref_id(const std::string &ref_label) const;
  //const std::map<id_t, std::string>& get_refs() const;
  const std::set<pt::id_t> &get_shared_samples(pt::id_t ref_id) const;
  pt::id_t ref_id_count() const;
  bool has_refs() const;

  // ---------
  // setter(s)
  // ---------
  void add_tip(pt::id_t v_id, pgt::v_end_e end);
  // returns the index (v_idx) of the added vertex
  pt::idx_t add_vertex(pt::id_t v_id, const std::string& label);
  // returns the index (e_idx) of the added edge
  pt::idx_t add_edge(pt::id_t v1_id, pgt::v_end_e v1_end, pt::id_t v2_id, pgt::v_end_e v2_end);
  pt::id_t add_ref(const std::string &label, char delim);
  void shrink_to_fit();

  // other
  void summary(bool print_tips) const;
  void print_dot(std::ostream& os) const;
};

typedef VariationGraph VG;

std::vector<VG *> componetize(const VG &g);

pst::Tree compute_spanning_tree(const VG &g);

/**
  * @brief Get the paths between the flubble start and end
 */
//void populate_walks(const VG &g, pvt::RoV &r, pt::idx_t max_steps);

} // namespace povu::bidirected
#endif
