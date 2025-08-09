#ifndef POVU_GENOMIC_TYPES_HPP
#define POVU_GENOMIC_TYPES_HPP

#include <format>
#include <map>
#include <memory>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

#include "./pvst.hpp"
#include "./graph.hpp"
#include "./core.hpp"

namespace povu::types::genomics {
namespace pgt = povu::types::graph;
namespace pt = povu::types;
namespace pc = povu::constants;
namespace pvst = povu::types::pvst;

/* ===========================================================
     types that involve the graph itself

      step, walk, RoV
   =========================================================== */

// TODO: move to graph types
typedef pgt::id_or_t step;
typedef std::vector<pgt::id_or_t> walk;

/**
 * a collection of walks within a region of variation from start to end
 */
class RoV {
std::vector<walk> walks_;
const pvst::VertexBase *pvst_vtx;

public:

// --------------
// constructor(s)
// --------------

RoV(const pvst::VertexBase *v) : walks_(), pvst_vtx(v) {}

// ---------
// getter(s)
// ---------

pt::idx_t walk_count() const { return this->walks_.size(); }
const pvst::VertexBase *get_pvst_vtx() const { return this->pvst_vtx; }
const std::vector<walk> &get_walks() const { return this->walks_; }
std::vector<walk> &get_walks_mut() { return this->walks_; }

// ---------
// setter(s)
// ---------

void set_walks(std::vector<walk> &&walks) { this->walks_ = walks; }

// --------
// other(s)
// --------

std::string as_str() const {
  return this->pvst_vtx->as_str();
}
};

/* ===========================================================
     types that involve the reference

      Allele Step, Allele Walk, Itinerary, Ref Walks
   =========================================================== */

enum class var_type_e {
  del, // deletion
  ins, // insertion
  sub, // substitution
  und // undetermined
};

// add operator << for var_type_e
inline std::ostream &operator<<(std::ostream &os, var_type_e vt) {
  switch (vt) {
    case var_type_e::del:
      os << "DEL";
      break;
    case var_type_e::ins:
      os << "INS";
      break;
    case var_type_e::sub:
      os << "SUB";
      break;
    case var_type_e::und:
      os << "UND";
      break;
  }
  return os;
}

// 1) A helper to turn enum → string_view:
constexpr std::string_view to_string_view(var_type_e vt) noexcept {
  switch (vt) {
  case var_type_e::del:
    return "DEL";
  case var_type_e::ins:
    return "INS";
  case var_type_e::sub:
    return "SUB";
  case var_type_e::und:
    return "UND";
  }
  // optional: handle out-of-range
  return "??";
}

/**
 * Allele or Ref Step
 */
class AS {
  // TODO: remove loop_no?
  // pt::id_t loop_no_; // the nth time that a ref is going through a flubble RoV
  pt::id_t v_id_;
  pt::idx_t step_idx_; // TODO: rename to locus_?
  pgt::or_e o_;

public:

  // --------------
  // constructor(s)
  // --------------
  AS(pt::id_t v_id, pgt::or_e o) :v_id_(v_id), step_idx_(pc::INVALID_IDX), o_(o) {}
  AS(pt::id_t v_id, pt::idx_t step_idx, pgt::or_e o ) :v_id_(v_id), step_idx_(step_idx), o_(o) {}

  // ---------
  // getter(s)
  // ---------
  pt::idx_t get_step_idx() const { return this->step_idx_; }
  pt::id_t get_v_id() const { return this->v_id_; }
  pgt::or_e get_o() const { return this->o_; }
  //pt::id_t get_loop_no() const { return this->loop_no_; }

  // ---------
  // setter(s)
  // ---------
  void set_step_idx(pt::idx_t step_idx) { this->step_idx_ = step_idx; }
  //void set_loop_no(pt::id_t loop_id) { this->loop_no_ = loop_id; }
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
  //const pvst::VertexBase *pvst_vtx; // TODO: remove
  //bool is_del_; // TODO: remove
  //var_type_e vt_;
  pt::idx_t walk_idx_;

  std::set<pt::id_t> ref_ids_; // set of ref ids that take this walk

public:

  // --------------
  // constructor(s)
  // --------------
  AW() : steps_() {}
  AW(pt::idx_t w_idx) : steps_(), walk_idx_(w_idx) {}
  AW(pt::id_t id, pgt::or_e o) : steps_(std::vector<AS>{AS{id, o}}) {}
  AW(AS s) : steps_(std::vector<AS>{s}) {}

  // ---------
  // getter(s)
  // ---------
  pt::idx_t step_count() const { return this->steps_.size(); }
  const std::vector<AS> &get_steps() const { return this->steps_; }
  std::vector<AS> &get_steps_mut() { return this->steps_; }
  const AS &get_step(pt::idx_t idx) const { return this->steps_[idx]; }
  std::string as_str() const {
    std::string s;
    for (const AS &step : this->steps_) {
      s += std::format("{}{}", step.get_o() == pgt::or_e::forward ? ">" : "<", step.get_v_id());
    }
    return s;
  }
  // bool is_del() const { return this->is_del_; }
  //var_type_e get_var_type() const { return this->vt_; }
  bool empty() const noexcept { return this->step_count() == 0; }
  pt::idx_t get_walk_idx() const { return this->walk_idx_; }
  const std::set<pt::id_t> &get_ref_ids() const { return this->ref_ids_; }


  // ---------
  // setter(s)
  // ---------
  void append_step(AS s) { this->steps_.emplace_back(s); }
  void add_ref_id(pt::id_t ref_id) { this->ref_ids_.insert(ref_id); }
  //void set_var_type(var_type_e vt) { this->vt_ = vt; }
  //void set_is_del(bool is_del) { this->is_del_ = is_del; }
};

// TODO: replace instances of at with aw

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

    this->it_[at_idx].append_step(s);
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
  std::map<pt::op_t<pt::id_t>, std::string> aln;

  const pvst::VertexBase *pvst_vtx;

  // is true when tangling exists.
  // tangling exists when a walk traverses an RoV more than once
  bool is_tangled_;

public:

  // --------------
  // constructor(s)
  // --------------

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
    return this->aln.at(pt::op_t<pt::id_t>{ref_id1, ref_id2});
  }

  bool has_aln(pt::id_t ref_id1, pt::id_t ref_id2) const {
    return this->aln.contains(pt::op_t<pt::id_t>{ref_id1, ref_id2});
  }

  const std::map<pt::op_t<pt::id_t>, std::string> &get_alns() const {
    return this->aln;
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
    this->aln[pt::op_t<pt::id_t>{ref_id1, ref_id2}] = aln;
  }

  void set_tangled(bool is_tangled) {
    this->is_tangled_ = is_tangled;
  }
};


class VcfRec {
  pt::id_t ref_id_; // chrom
  pt::idx_t pos; // 1-based step idx
  std::string id; // start and end of a RoV e.g >1>4
  AW ref_at; //
  std::vector<AW> alt_ats;
  // std::string qual;
  // std::string filter;
  // std::string info;
  std::string format;
  pt::idx_t height_; // height of the pvst node in the tree
  var_type_e var_type_; // type of the variant, e.g. del, ins, sub, und
  bool is_tangled_ = false; // is true when tangling exists, i.e. when a walk traverses an RoV more than once

public:

  // --------------
  // constructor(s)
  // --------------

  VcfRec(pt::id_t ref_id, pt::idx_t pos, std::string id, AW ref_at,
         std::vector<AW> alt_ats, pt::idx_t height, var_type_e variant_type, bool is_tangled)
    : ref_id_(ref_id), pos(pos), id(id), ref_at(ref_at), alt_ats(alt_ats),
      height_(height), var_type_(variant_type), is_tangled_(is_tangled) {}

  // ---------
  // getter(s)
  // ---------
  pt::id_t get_ref_id() const { return this->ref_id_; }
  pt::idx_t get_pos() const { return this->pos; }
  std::string get_id() const { return this->id; }
  const AW &get_ref_at() const { return this->ref_at; }
  const std::vector<AW> &get_alt_ats() const { return this->alt_ats; }
  std::vector<AW> &get_alt_ats_mut() { return this->alt_ats; }
  pt::idx_t get_height() const { return this->height_; }
  var_type_e get_var_type() const { return this->var_type_; }
  bool is_tangled() const { return this->is_tangled_; }

  std::set<pt::id_t> get_ref_ids() const {
    std::set<pt::id_t> ref_ids;
    ref_ids.insert(this->ref_id_);
    for (const auto &aw : this->alt_ats) {
      for (const auto &ref_id : aw.get_ref_ids()) {
        ref_ids.insert(ref_id);
      }
    }
    return ref_ids;
  }


  // ---------
  // setter(s)
  // ---------
  pt::idx_t append_alt_at(AW &&alt_at) {
    this->alt_ats.emplace_back(std::move(alt_at));
    return this->alt_ats.size() - 1;
  }
};


struct genotype_data_t {
  std::map<pt::id_t, pt::idx_t> ref_id_to_col_idx;

  // the size is the number of columns in the genotype data
  // string contains the sample name or label
  std::vector<std::string> genotype_cols;
};

// TODO: [c] find a better name
class VcfRecIdx {
  std::map<pt::idx_t, std::vector<VcfRec>> vcf_recs_;
  genotype_data_t gd_;

public:

  // --------------
  // constructor(s)
  // --------------

  VcfRecIdx() : vcf_recs_() {}

  // ---------
  // getter(s)
  // ---------
  const std::map<pt::idx_t, std::vector<VcfRec>> &get_recs() const {
    return this->vcf_recs_;
  }

  const genotype_data_t &get_genotype_data() const {
    return this->gd_;
  }

  // ---------
  // setter(s)
  // ---------

  void add_rec(pt::idx_t ref_idx, VcfRec &&vcf_rec) {
    if (this->vcf_recs_.find(ref_idx) == this->vcf_recs_.end()) {
      this->vcf_recs_[ref_idx] = std::vector<VcfRec>{vcf_rec};
    } else {
      this->vcf_recs_[ref_idx].emplace_back(vcf_rec);
    }
  }

  void set_genotype_data(genotype_data_t &&gd) {
    this->gd_ = std::move(gd);
  }

};


/* ===========================================================
     enums

      aln_level_e
   =========================================================== */

/**
 * used in untangling, the level of alignment
 */
enum class aln_level_e {
  step,
  at // allele traversal
};

} // namespace povu::types::genomics
#endif // POVU_GENOMIC_TYPES_HPP
