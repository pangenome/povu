#ifndef POVU_GENOMIC_TYPES_HPP
#define POVU_GENOMIC_TYPES_HPP

#include <map>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "./types.hpp"


namespace povu::types::genomics {
namespace pgt = povu::graph_types;
namespace pt = povu::types;
namespace pc = povu::constants;

class Step {
  pt::id_t loop_no_; // the nth time that a ref is going through a flubble RoV
  pt::id_t v_id_;
  pt::idx_t step_idx_;
  pgt::or_e o_;

public:
/* constructor */
Step(pt::id_t v_id, pgt::or_e o)
  : loop_no_(pc::INVALID_ID), v_id_(v_id), step_idx_(pc::INVALID_IDX), o_(o) {}

/*getters*/
pt::idx_t get_step_idx() const { return this->step_idx_; }
pt::id_t get_v_id() const { return this->v_id_; }
pgt::or_e get_o() const { return this->o_; }
pt::id_t get_loop_no() const { return this->loop_no_; }

/*setters*/
void set_step_idx(pt::idx_t step_idx) { this->step_idx_ = step_idx; }
void set_loop_no(pt::id_t loop_id) { this->loop_no_ = loop_id; }
};


/* an uninterrupted ordered sequence of steps bound by the start end of a RoV */
class Walk {
std::vector<Step> steps_;

public:
/* constructor */
Walk() : steps_() {}
Walk(pt::id_t id, pgt::or_e o) : steps_(std::vector<Step>{Step{id, o}}) {}

/*getters*/
pt::idx_t step_count() const { return this->steps_.size(); }
const std::vector<Step> &get_steps() const { return this->steps_; }
std::vector<Step> &get_steps_mut() { return this->steps_; }
const Step &get_step(pt::idx_t idx) const { return this->steps_[idx]; }
std::string as_str() const {
  std::string s;
    for (const Step &step : this->steps_) {

      s += std::format("{}{}", step.get_o() == pgt::or_e::forward ? ">" : "<",
                       step.get_v_id());
    }
    return s;
  }

/*setters*/
void append_step(Step s) { this->steps_.emplace_back(s); }
};

typedef Walk StepSeq;

/*allele traversal---a walk taken by a reference*/
typedef Walk AT;

/* region of variation */
class RoV {
std::vector<Walk> walks_;
pgt::flubble_t fl_;

public:
RoV(pgt::flubble_t fl) : walks_(), fl_(fl) {}

/*getters*/
pt::idx_t walk_count() const { return this->walks_.size(); }
const pgt::id_or_t &get_entry() const { return this->fl_.start_; }
const pgt::id_or_t &get_exit() const { return this->fl_.end_; }
pgt::flubble_t get_flb() const { return this->fl_; }
const std::vector<Walk> &get_walks() const { return this->walks_; }
std::vector<Walk> &get_walks_mut() { return this->walks_; }

/*setters*/
void set_walks(std::vector<Walk> &&walks) { this->walks_ = walks; }

/* other */
std::string as_str() const {
  return this->fl_.start_.as_str() + this->fl_.end_.as_str();
}
};

/*
 * an interrupted
 * sequence of looped walks in a RoV for a given ref
*/
class It {
  std::vector<AT> it_;
  pt::idx_t len;

public:
  /* constructors */
  It() : it_(), len(0) {}
  It(AT &&w) : it_(std::vector<AT>{w}), len(w.step_count()) {}

  /* getters */
  pt::idx_t walk_count() const { return this->it_.size(); }
  pt::idx_t step_count() const { return this->len; }

  const std::vector<AT> &get_walks() const { return this->it_; }
  std::vector<AT> &get_walks_mut() { return this->it_; }

  AT &get_walk_mut(pt::idx_t w_idx) { return this->it_[w_idx]; }

  const Step &get_step(pt::idx_t idx) const {
    pt::idx_t i = 0;
    for (const AT &w : this->it_) {
      if (i + w.step_count() > idx) {
        return w.get_steps()[idx - i];
      }
      i += w.step_count();
    }
    throw std::out_of_range("step index out of range");
  }

  const AT &get_walk_by_step_idx(pt::idx_t idx) const {
    pt::idx_t i = 0;
    for (const AT &w : this->it_) {
      if (i + w.step_count() > idx) {
        return w;
      }
      i += w.step_count();
    }
    throw std::out_of_range("step index out of range");
  }

  /* setters */
  void add_walk(AT &&w) {
    this->it_.emplace_back(w);
    this->len += w.step_count();
  }
};

typedef It Itn;

/*  map of ref_id to the walk of the ref in a RoV */
class RefWalks {
  // up = unordered pair
  typedef std::pair<pt::id_t, pt::id_t> up_t;

  // use Walk instead of vector<Step>?
  // map of ref_id to the walk of the ref in a RoV
  // a ref can have multiple walks in a RoV
  // the order of walks is not assured
  std::map<pt::id_t, Itn> ref_walks_;

  //alignment between two refs
  std::map<up_t, std::string> aln;

  pgt::flubble_t fl_;

  bool is_tangled_;

  /* private methods */
  // returns an unordered pair
  up_t to_up (pt::id_t a, pt::id_t b) const {
    return {std::min(a, b), std::max(a, b)};
  };

public:
  /* constructor */
  RefWalks(pgt::flubble_t fl) : ref_walks_(), fl_(fl) {
    is_tangled_ = false;
  }

  /*getters*/
  pt::idx_t ref_count() const { return this->ref_walks_.size(); }
  const pgt::flubble_t &get_flb() const { return this->fl_; }

  std::set<pt::id_t> get_ref_ids() const {
    std::set<pt::id_t> ref_ids;
    for (const auto &p : this->ref_walks_) {
      ref_ids.insert(p.first);
    }
    return ref_ids;
  }

  const Itn &get_itn(pt::id_t ref_id) const {
    return this->ref_walks_.at(ref_id);
  }

  const std::map<pt::id_t, It> &get_ref_walks() const {
    return this->ref_walks_;
  }

  bool has_aln(pt::id_t ref_id1, pt::id_t ref_id2) const {
    return this->aln.find(to_up(ref_id1, ref_id2)) != this->aln.end();
  }

  const std::string &get_aln(pt::id_t ref_id1, pt::id_t ref_id2) const {
    return this->aln.at(to_up(ref_id1, ref_id2));
  }

  const std::map<up_t, std::string> &get_alns() const {
    return this->aln;
  }

  bool is_tangled() const { return this->is_tangled_; }

  /*setters*/
  // returns true to mean that tangling exists. a walk traverses an RoV more than once
  void add_itn(id_t ref_id, Itn &&itn) {
    if (this->ref_walks_.find(ref_id) == this->ref_walks_.end()) {
      this->ref_walks_[ref_id] = std::move(itn);
    } else {
      Itn &itn_ = this->ref_walks_[ref_id];
      for (auto &w : itn.get_walks_mut()) {
        itn_.add_walk(std::move(w));
      }

      is_tangled_ = true;
    }
  }


  void sort_by_step_idx() {
    // for each ref, sort the walks by step_idx of the first step
    for (auto &[_, itn] : this->ref_walks_) {
      std::sort(itn.get_walks_mut().begin(), itn.get_walks_mut().end(),
                [](const Walk &a, const Walk &b) {
                  return a.get_steps().front().get_step_idx() <
                         b.get_steps().front().get_step_idx();
                });
    }

    // set the loop number for each step
    // steps within the same walk share the same loop number
    for (auto &[_, itn] : this->ref_walks_) {
      pt::id_t loop_no = 0;
      for (Walk &w : itn.get_walks_mut()) {
        for (Step &s : w.get_steps_mut()) {
          s.set_loop_no(loop_no);
        }
        loop_no++;
      }
    }
  }

  void add_aln(pt::id_t ref_id1, pt::id_t ref_id2 , std::string &&aln) {
    this->aln[to_up(ref_id1, ref_id2)] = aln;
  }
};

class VcfRec {
  pt::id_t ref; // chrom
  pt::idx_t pos; // 1-based step idx
  std::string id; // start and end of a RoV e.g >1>4
  AT ref_at; //
  std::vector<AT> alt_ats;
  // std::string qual;
  // std::string filter;
  // std::string info;
  std::string format;

public:
  /* constructor */
  VcfRec(pt::id_t ref, pt::idx_t pos, std::string id, Walk ref_at,
         std::vector<Walk> alt_ats, std::string format)
      : ref(ref), pos(pos), id(id), ref_at(ref_at), alt_ats(alt_ats), format(format) {}

  /*getters*/
  pt::idx_t get_pos() const { return this->pos; }
  std::string get_id() const { return this->id; }
  const AT &get_ref_at() const { return this->ref_at; }
  const std::vector<AT> &get_alt_ats() const { return this->alt_ats; }

  /*setters*/
};

// TODO: find a better name
class VcfRecIdx {
  std::map<pt::idx_t, std::vector<VcfRec>> vcf_recs_;

  public:
  /* constructor */
  VcfRecIdx() : vcf_recs_() {}

  /*getters*/
  const std::map<pt::idx_t, std::vector<VcfRec>> &get_recs() const {
    return this->vcf_recs_;
  }

  /*setters*/
  void add_rec(pt::idx_t ref_idx, VcfRec &&vcf_rec) {
    if (this->vcf_recs_.find(ref_idx) == this->vcf_recs_.end()) {
      this->vcf_recs_[ref_idx] = std::vector<VcfRec>{vcf_rec};
    } else {
      this->vcf_recs_[ref_idx].emplace_back(vcf_rec);
    }
  }

};
/*
  enums
  -----
*/
enum class aln_level_e {
  step,
  at // allele traversal
};

} // namespace povu::types::variation

#endif // POVU_GENOMIC_TYPES_HPP
