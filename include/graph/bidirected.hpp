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
#include <sys/types.h>
#include <tuple>
#include <unordered_set>
#include <vector>

#include "../../app/cli/app.hpp" // for core::config TODO: find a proper place for the app config
#include "../common/compat.hpp"
#include "../common/types/constants.hpp"
#include "../common/types/graph.hpp"
#include "../common/utils.hpp"
#include "../common/log.hpp"


namespace povu::bidirected {
inline constexpr std::string_view MODULE = "povu::bidirected";

namespace pu = povu::utils;
namespace pc = povu::constants;
using namespace povu::types::graph;
namespace pgt = povu::types::graph;

// TODO: come up with clear definitions for:
// haplotype, ref, contig, path, etc.
class RefInfo {
  pt::id_t ref_id_;
  pgt::or_e strand_; // TODO [c] chose between names orientation and strand
  pt::idx_t locus_;

public:
  // --------------
  // constructor(s)
  // --------------
  RefInfo(pt::id_t ref_id, pgt::or_e strand, pt::idx_t locus);

  // ---------
  // getter(s)
  // ---------
  pt::id_t get_ref_id() const;
  pgt::or_e get_strand() const;
  pt::idx_t get_locus() const;
};


class VtxRefIdx {
  // the set of ref ids that are present in the vertex
  std::set<pt::idx_t> refs_;
  std::set<pt::idx_t> loci_;

  // ref id to the index of the ref in the RefInfo vector
  std::map<pt::id_t, std::set<pt::idx_t>> ref_id_to_ref_info_idx_;
  // key is the ref id and value is a map of locus to the index of the ref in
  // the RefInfo vector
  std::map<pt::id_t, std::map<pt::idx_t, pt::idx_t>> ref_id_to_locus_idx_;

  // map of ref id to the set of loci where the ref is present in the vertex
  std::map<pt::id_t, std::set<pt::idx_t>> ref_id_to_loci_;

  // --------------
  // constructor(s)
  // --------------

  VtxRefIdx(){}

public :
  // --------------
  // factory method(s)
  // --------------

  static VtxRefIdx from_ref_info(std::vector<RefInfo> &v_ref_data) {
    VtxRefIdx vr_idx;
    for (pt::idx_t i = 0; i < v_ref_data.size(); ++i) {
      const RefInfo &ref = v_ref_data[i];
      pt::id_t ref_id = ref.get_ref_id();
      pt::idx_t locus = ref.get_locus();

      vr_idx.ref_id_to_ref_info_idx_[ref_id].insert(i);
      vr_idx.ref_id_to_locus_idx_[ref_id][locus] = i;
      vr_idx.ref_id_to_loci_[ref_id].insert(locus);
      vr_idx.refs_.insert(ref_id);
      vr_idx.loci_.insert(locus);
    }
    return vr_idx;
  }

  // --------------
  // getter(s)
  // --------------
  bool has_locus(pt::id_t ref_id, pt::idx_t qry_locus) const {
    if (!pv_cmp::contains(this->refs_, ref_id)) {
      return false;
    }
    const std::set<pt::idx_t> &loci_set = this->ref_id_to_loci_.at(ref_id);
    return pv_cmp::contains(loci_set, qry_locus);
  }

  pt::idx_t get_ref_data_idx(pt::id_t ref_id, pt::idx_t qry_locus) const {
    if (!pv_cmp::contains(this->ref_id_to_locus_idx_, ref_id)) {
      return pc::INVALID_IDX;
    }
    if (!pv_cmp::contains(this->ref_id_to_locus_idx_.at(ref_id), qry_locus)) {
      return pc::INVALID_IDX;
    }
    return this->ref_id_to_locus_idx_.at(ref_id).at(qry_locus);
  }

  const std::set<pt::idx_t> &get_ref_ids() const { return this->refs_; }

  const std::set<pt::idx_t> &get_ref_loci(pt::id_t ref_id) const {
    if (!pv_cmp::contains(this->ref_id_to_loci_, ref_id)) {
      ERR("ref id {} not found in vertex\n", ref_id);
      exit(1);
    }
    if (this->ref_id_to_loci_.at(ref_id).empty()) {
      ERR("no loci for ref {} in vertex\n", ref_id);
      exit(1);
    }
    return this->ref_id_to_loci_.at(ref_id);
  }

  // works with 0 because we assume indexes are 1 indexed
  // the min must be > threshold
  std::optional<pt::idx_t> get_min_locus(pt::id_t ref_id, std::optional<pt::idx_t> opt_start_after) const {
    if (!pv_cmp::contains(this->refs_, ref_id)) {
      return std::nullopt;
    }
    if (this->ref_id_to_loci_.at(ref_id).empty()) {
      ERR("no loci for ref {} in vertex\n", ref_id);
      exit(1);
    }

    pt::idx_t threshold = opt_start_after.value_or(0);
    for (pt::idx_t locus : this->ref_id_to_loci_.at(ref_id)) {
      if (locus > threshold) {
        return locus;
      }
    }
    return std::nullopt;
  }

  pt::idx_t loop_count(pt::id_t ref_id) const {
    if (!pv_cmp::contains(this->refs_, ref_id)) {
      return 0;
    }
    return this->ref_id_to_loci_.at(ref_id).size();
  }
};


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
  std::vector<RefInfo> refs_;

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
  const std::vector<RefInfo>& get_refs() const;

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


  std::map<pt::id_t, VtxRefIdx> v_ref_idx_; // v_id to VtxRefIdx
  bool has_refs_;
  pgt::Refs refs_;
  std::set<pgt::side_n_id_t> tips_; // the set of side and id of the tips



public:
  // --------------
  // constructor(s)
  // --------------
  VariationGraph(pt::idx_t vtx_count, pt::idx_t edge_count, bool inc_refs);

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
  const std::string &get_ref_label(pt::id_t ref_id) const;
  const pgt::Ref &get_ref_by_id(pt::id_t ref_id) const;
  pgt::Ref &get_ref_by_id_mut(pt::id_t ref_id);
  pt::id_t get_ref_id(const std::string &ref_label) const;
  //const std::map<id_t, std::string>& get_refs() const;
  const std::set<pt::id_t> &get_shared_samples(pt::id_t ref_id) const;
  pt::id_t ref_id_count() const;
  bool has_refs() const;
  const VtxRefIdx &get_vtx_ref_idx(pt::id_t v_id) const;

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

  void gen_vtx_ref_idx(pt::id_t v_id);

  // other
  void summary(bool print_tips) const;
  void print_dot(std::ostream& os) const;
};

typedef VariationGraph VG;
} // namespace povu::bidirected

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace bd = povu::bidirected;

#endif
