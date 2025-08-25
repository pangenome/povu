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

#include "../common/types/graph.hpp"
#include "../graph/pvst.hpp"
#include "../common/log.hpp"
#include "../graph/bidirected.hpp"

namespace povu::genomics::allele {
inline constexpr std::string_view MODULE = "povu::genomics::allele";

namespace pgt = povu::types::graph;
namespace pvst = povu::pvst;

/**
 * Allele or Ref Step
 */
class AS {
  pt::id_t v_id_;
  pt::idx_t step_idx_; // TODO: rename to locus_?
  pgt::or_e o_;

  // --------------
  // constructor(s)
  // --------------
  AS(pt::id_t v_id, pt::idx_t step_idx, pgt::or_e o)
    : v_id_(v_id), step_idx_(step_idx), o_(o) {}

public:

  // -----------------
  // factory method(s)
  // -----------------
  static AS given_ref_info(pt::id_t v_id, const bd::RefInfo &ref_info) {
    return AS(v_id, ref_info.get_locus(), ref_info.get_strand());
  }

  // ---------
  // getter(s)
  // ---------
  pt::idx_t get_step_idx() const { return this->step_idx_; }
  pt::id_t get_v_id() const { return this->v_id_; }
  pgt::or_e get_o() const { return this->o_; }

  // ---------
  // setter(s)
  // ---------
  void set_step_idx(pt::idx_t step_idx) { this->step_idx_ = step_idx; }
};

// implement == operator for AS
inline bool operator==(const AS &lhs, const AS &rhs) {
  return lhs.get_v_id() == rhs.get_v_id() &&
         lhs.get_step_idx() == rhs.get_step_idx() &&
         lhs.get_o() == rhs.get_o();
}

/**
 * Ref walk or Allele Walk
 * allele walk---a walk taken by a reference
 * an uninterrupted ordered sequence of steps bound by the start end of a RoV
 */
class AW {
  std::vector<AS> steps_;
  pt::idx_t walk_idx_;
  std::set<pt::id_t> ref_ids_; // set of ref ids that take this walk

public:

  // --------------
  // constructor(s)
  // --------------
  //AW() : steps_() {}
    AW(pt::idx_t w_idx) : steps_(), walk_idx_(w_idx) {}
  //AW(AS s) : steps_(std::vector<AS>{s}) {}


  // ---------
  // getter(s)
  // ---------
  pt::idx_t step_count() const { return this->steps_.size(); }
  const std::vector<AS> &get_steps() const { return this->steps_; }
  std::vector<AS> &get_steps_mut() { return this->steps_; }
  const AS &get_step(pt::idx_t idx) const { return this->steps_[idx]; }
  // TODO: use fmt format for user defined types
  // https://fmt.dev/11.1/api/#udt
  std::string as_str() const {
    std::string res;
    for (const AS &s : this->steps_) {
      // TODO: [C] [PERF] replace with format_to?
      // TODO: use some kind of string formatting for pgt::orientation
      res += pv_cmp::format("{}{}", s.get_o() == pgt::or_e::forward ? ">" : "<", s.get_v_id());
    }
    return res;
  }
  bool empty() const noexcept { return this->step_count() == 0; }
  pt::idx_t get_walk_idx() const { return this->walk_idx_; }
  const std::set<pt::id_t> &get_ref_ids() const { return this->ref_ids_; }


  // ---------
  // setter(s)
  // ---------
  void append_step(AS &&s) { this->steps_.emplace_back(s); }
  void add_ref_id(pt::id_t ref_id) { this->ref_ids_.insert(ref_id); }
};


/**
 * Ref Itinerary or just Itinerary
 * an interrupted
 * sequence of looped walks in a RoV for a given ref
 * useful for repeats
 */
class Itn {
  std::vector<AW> it_;
  pt::idx_t len; // total step count in the itinerary

public:

  // --------------
  // constructor(s)
  // --------------
  Itn() : it_(), len(0) {}

  // ---------
  // getter(s)
  // ---------
  pt::idx_t at_count() const { return this->it_.size(); }
  pt::idx_t step_count() const { return this->len; }

  const std::vector<AW> &get_ats() const { return this->it_; }
  std::vector<AW> &get_ats_mut() { return this->it_; }

  const AW &get_at(pt::idx_t at_idx) const { return this->it_[at_idx]; }
  AW &get_at_mut(pt::idx_t at_idx) { return this->it_[at_idx]; }

  const AS &get_step(pt::idx_t idx) const {
    pt::idx_t i = 0;
    for (const AW &w : this->it_) {
      if (i + w.step_count() > idx) {
        return w.get_steps()[idx - i];
      }
      i += w.step_count();
    }
    throw std::out_of_range("step index out of range");
  }

  // ---------
  // setter(s)
  // ---------
  void append_at(AW &&w) {
    this->it_.emplace_back(w);
    this->len += w.step_count();
  }

  void append_step(pt::idx_t at_idx, AS s) {
    if (at_idx >= this->at_count()) {
      throw std::out_of_range("at index out of range");
    }

    this->it_[at_idx].append_step(std::move(s));
    this->len++;
  }

  void remove_at(pt::idx_t at_idx) {
    if (at_idx >= this->at_count()) {
      throw std::out_of_range("at index out of range");
    }

    this->len -= this->it_[at_idx].step_count();
    this->it_.erase(this->it_.begin() + at_idx);
  }

  // If you delete the highest indices first, the lower-indexed ones never
  // shift:
  void remove_aws(const std::set<pt::idx_t> &to_remove) {
    // walk the set in reverse (largest index → smallest)
    for (auto it = to_remove.rbegin(); it != to_remove.rend(); ++it) {
      pt::idx_t aw_idx = *it;
      if (aw_idx >= this->at_count())
        throw std::out_of_range("…");

      // now it's safe to use aw_idx directly:
      this->len -= this->it_[aw_idx].step_count();
      this->it_.erase(this->it_.begin() + aw_idx);
    }
  }
};

/**
 * Exp, short for expedition: a journey undertaken by a group of people with
 * a particular purpose, especially that of exploration, research, or war.
 *
 * A collection of itineraries for each reference in a region of variation
 * map of ref_id to the itn of the ref in a RoV
 */
class Exp {
  // map of ref_id to the itinerary (set of walks) of the ref in a RoV
  // when tangled, a ref can have multiple walks in a RoV
  std::map<pt::id_t, Itn> ref_itns_;

  // alignment between two refs
  std::map<pt::op_t<pt::id_t>, std::string> aln_;

  const pvst::VertexBase *pvst_vtx;

  // is true when tangling exists.
  // tangling exists when a walk traverses an RoV more than once
  bool is_tangled_;

public:

  // --------------
  // constructor(s)
  // --------------
  Exp() : ref_itns_(), pvst_vtx(nullptr), is_tangled_(false) {}
  Exp(const pvst::VertexBase *v) : ref_itns_(), pvst_vtx(v), is_tangled_(false) {}

  // ---------
  // getter(s)
  // ---------

  pt::idx_t ref_count() const { return this->ref_itns_.size(); }

  const pvst::VertexBase* get_pvst_vtx_const_ptr() const {
    return this->pvst_vtx;
  }

  std::string id() const {
    return this->pvst_vtx->as_str();
  }

  std::set<pt::id_t> get_ref_ids() const {
    std::set<pt::id_t> ref_ids;
    for (const auto &p : this->ref_itns_) {
      ref_ids.insert(p.first);
    }
    return ref_ids;
  }

  const Itn &get_itn(pt::id_t ref_id) const {
    return this->ref_itns_.at(ref_id);
  }

  Itn &get_itn_mut(pt::id_t ref_id) {
    return this->ref_itns_.at(ref_id);
  }

  const std::map<pt::id_t, Itn> &get_ref_itns() const {
    return this->ref_itns_;
  }

  std::map<pt::id_t, Itn> &get_ref_itns_mut() {
    return this->ref_itns_;
  }

  const std::string &get_aln(pt::id_t ref_id1, pt::id_t ref_id2) const {
    return this->aln_.at(pt::op_t<pt::id_t>{ref_id1, ref_id2});
  }

  bool has_aln(pt::id_t ref_id1, pt::id_t ref_id2) const {
    return pv_cmp::contains(this->aln_, pt::op_t<pt::id_t>{ref_id1, ref_id2});
  }

  const std::map<pt::op_t<pt::id_t>, std::string> &get_alns() const {
    return this->aln_;
  }

  bool is_tangled() const { return this->is_tangled_; }

  // ---------
  // setter(s)
  // ---------

  void add_itn(id_t ref_id, Itn &&itn) {
    if (this->ref_itns_.find(ref_id) == this->ref_itns_.end()) {
      this->ref_itns_[ref_id] = std::move(itn);
    }
    else {
      Itn &itn_ = this->ref_itns_[ref_id];
      for (auto &w : itn.get_ats_mut()) {
        itn_.append_at(std::move(w));
      }
    }
  }

  void replace_itn(pt::id_t ref_id, Itn &&itn) {
    this->ref_itns_[ref_id] = std::move(itn);
  }

  void add_aln(pt::id_t ref_id1, pt::id_t ref_id2 , std::string &&aln) {
    this->aln_[pt::op_t<pt::id_t>{ref_id1, ref_id2}] = aln;
  }

  void set_tangled(bool is_tangled) {
    this->is_tangled_ = is_tangled;
  }
};


class WalkRefIdx {
  std::set<pt::idx_t> refs_;

  // --------------
  // constructor(s)
  // --------------
  WalkRefIdx() {}

public:
  // --------------
  // factory method(s)
  // --------------
  static WalkRefIdx from_walk(const bd::VG &g, const pgt::walk_t &w) {
    WalkRefIdx wr_idx;
    for (const auto &[v_id, _] : w) {
      const bd::VtxRefIdx &vr_idx = g.get_vtx_ref_idx(v_id);
      wr_idx.refs_.insert(vr_idx.get_ref_ids().begin(), vr_idx.get_ref_ids().end());
    }

    return wr_idx;
  }

  // --------------
  // getter(s)
  // --------------
  const std::set<pt::idx_t> &get_ref_ids() const { return this->refs_; }
};

void comp_itineraries(const bd::VG &g, const std::vector<pgt::walk_t> &walks,
                      std::map<pt::id_t, Itn> &ref_map);

} // namespace povu::genomics::alele

#endif // POVU_GENOMICS_ALLELE_HPP
