#ifndef BIDIRECTED_HPP
#define BIDIRECTED_HPP

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <iterator>
#include <map>
#include <optional>
#include <ostream>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <string_view>
#include <sys/types.h>
#include <tuple>
#include <unordered_set>
#include <vector>


#include "../common/compat.hpp"
#include "../common/constants.hpp"
#include "./types.hpp"
#include "../common/utils.hpp"
#include "../common/log.hpp"
#include "../refs/refs.hpp"
#include "./ref.hpp"

namespace povu::bidirected {
inline constexpr std::string_view MODULE = "povu::bidirected";

namespace pu = povu::utils;
namespace pc = povu::constants;
using namespace povu::types::graph;
namespace pgt = povu::types::graph;


// undirected edge
// stores the index of the vertex in the graph not the id
class Edge {
  pt::idx_t v1_idx_;
  pgt::v_end_e v1_end_;
  pt::idx_t v2_idx_;
  pgt::v_end_e v2_end_;

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
  pt::id_t v_id_; // this is the sequence name in the GFA file. Maybe Should support strings.
  std::string label_; // or sequence

  // indexes to the edge vector in Graph
  std::set<pt::idx_t> e_l;
  std::set<pt::idx_t> e_r;

  // references (also colors)
  pgr::VtxRefData refs_;

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
  const pgr::VtxRefData& get_refs() const;

  // ---------
  // setter(s)
  // ---------
  void add_edge_l(pt::idx_t e_idx);
  void add_edge_r(pt::idx_t e_idx);
  void add_ref(pt::id_t ref_id, pgt::or_e strand, pt::idx_t locus);
};


class VariationGraph {
  std::vector<Vertex> vertices;
  std::vector<Edge> edges;
  pu::TwoWayMap<std::size_t, std::size_t> v_id_to_idx_; // TODO: reserve size

  pr::Refs refs_ = pr::Refs(0); // has no refs by default
  std::set<pgt::side_n_id_t> tips_; // the set of side and id of the tips

public:
  // --------------
  // constructor(s)
  // --------------
  VariationGraph(pt::idx_t vtx_count, pt::idx_t edge_count, pt::idx_t ref_count);

  // -----------------
  // factory method(s)
  // -----------------
  // return a vector of connected components as VG objects
  static std::vector<VariationGraph *> componetize(const VariationGraph &g);

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
  const std::string &get_sample_name(pt::id_t ref_id) const;
  const pr::Ref &get_ref_by_id(pt::id_t ref_id) const;
  pr::Ref &get_ref_by_id_mut(pt::id_t ref_id);
  std::optional<pt::id_t> get_ref_id(std::string_view ref_tag) const;
  pt::id_t ref_count() const;
  bool has_refs() const;
  /**
   * if PanSN the prefix will be in sample
   */
  std::set<pt::id_t> get_shared_samples(pt::id_t ref_id) const;
  // sometimes the sample name is referred to as a prefix
  std::set<pt::id_t> get_refs_in_sample(std::string_view sample_name) const;



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


  // -----
  // other
  // -----
  void summary(bool print_tips) const;
  void print_dot(std::ostream& os) const;
};

typedef VariationGraph VG;
} // namespace povu::bidirected

// NOLINTNEXTLINE(misc-unused-alias-decls)

namespace bd = povu::bidirected;

#endif
