#ifndef POVU_GENOMICS_ALLELE_HPP
#define POVU_GENOMICS_ALLELE_HPP

#include <algorithm>
#include <deque>
#include <functional> // for std::reference_wrapper
#include <future>
#include <optional>


#include <set>
#include <stack>
#include <string>
#include <sys/types.h>
#include <utility>
#include <vector>

#include "../common/compat.hpp"
#include "../common/types/core.hpp"
#include "../common/types/genomics.hpp"
#include "../common/types/graph.hpp"
#include "../common/types/pvst.hpp"
#include "../common/log.hpp"
#include "../graph/bidirected.hpp"

namespace povu::genomics::allele {
inline constexpr std::string_view MODULE = "povu::genomics::allele";

namespace pvt = povu::types::genomics;


// the ref visits of a single vertex. Unlike the ref visits in a walk, this
// contains the steps that the ref takes in the vertex, sorted by step index
typedef std::vector<pvt::AS> VtxRefVisits;


class VtxRefIdx {
  // ref id and start loci for this ref in the vertex
  //std::set<std::pair<pt::idx_t, pt::idx_t>> d_;

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

  static VtxRefIdx from_ref_info(std::vector<bd::RefInfo> &v_ref_data) {
    VtxRefIdx vr_idx;
    for (pt::idx_t i = 0; i < v_ref_data.size(); ++i) {
      const bd::RefInfo &ref = v_ref_data[i];
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

  std::optional<std::reference_wrapper<const std::set<pt::idx_t>>>
  get_ref_loci(pt::id_t ref_id) const {
    if (!pv_cmp::contains(this->ref_id_to_loci_, ref_id)) {
      return std::nullopt;
    }
    if (this->ref_id_to_loci_.at(ref_id).empty()) {
      ERR("no loci for ref {} in vertex\n", ref_id);
      exit(1);
    }
    return this->ref_id_to_loci_.at(ref_id);
  }

  };

class WalkRefIdx {
  // v_id to VtxRefIdx
  std::map<pt::id_t, VtxRefIdx> walk_ref_meta_;
  std::set<pt::idx_t> refs_;

  // --------------
  // constructor(s)
  // --------------
  WalkRefIdx() {}

public:
  // --------------
  // factory method(s)
  // --------------
  static WalkRefIdx from_walk(const bd::VG &g, const pvt::walk_t &w) {

    // handle dups
    std::set<pvt::step_t> seen;
    auto is_dup = [&](const pvt::step_t &s) -> bool {
      if (pv_cmp::contains(seen, s)) {
        WARN("step {} seen twice in walk \n", s.as_str());
        return true;
      }
      return false;
    };

    WalkRefIdx wr_idx;
    for (const pvt::step_t &s : w) {
      if (is_dup(s)) {
        continue;
      }

      pt::id_t v_id = s.v_id;
      const bd::Vertex &v = g.get_vertex_by_id(v_id);
      std::vector<bd::RefInfo> v_ref_data = v.get_refs();

      VtxRefIdx vr_idx = VtxRefIdx::from_ref_info(v_ref_data);
      wr_idx.refs_.insert(vr_idx.get_ref_ids().begin(), vr_idx.get_ref_ids().end());
      wr_idx.walk_ref_meta_.emplace(v_id, std::move(vr_idx));
    }

    return wr_idx;
  }

  // --------------
  // getter(s)
  // --------------
  const std::set<pt::idx_t> &get_ref_ids() const { return this->refs_; }


  const VtxRefIdx &get_vtx_ref_idx(pt::id_t v_id) const {
    return this->walk_ref_meta_.at(v_id);
  }

  pt::idx_t loop_count(pt::id_t v_id, pt::id_t ref_id) const {
    if (!pv_cmp::contains(this->refs_, ref_id)) {
      return 0;
    }
    const VtxRefIdx &vri = this->walk_ref_meta_.at(v_id);
    if (auto opt_loci = vri.get_ref_loci(ref_id); opt_loci) {
      return opt_loci->get().size(); // the set is expected to be non-empty
    }

    return 0;
  }

  std::optional<pt::idx_t> get_min_locus(pt::idx_t v_id, pt::idx_t ref_id) const {
    if (!pv_cmp::contains(this->refs_, ref_id)) {
      return std::nullopt;
    }
    const VtxRefIdx &vri = this->walk_ref_meta_.at(v_id);
    if (auto opt_loci = vri.get_ref_loci(ref_id); opt_loci) {
      return *opt_loci->get().begin(); // the set is expected to be non-empty
    }

    return std::nullopt;
  }

  // --------------
  // setters
  // --------------


};




void comp_itineraries(const bd::VG &g, const pvt::walk_t &w, pt::idx_t w_idx, pvt::Exp &rw);

} // namespace povu::genomics::alele

#endif // POVU_GENOMICS_ALLELE_HPP
