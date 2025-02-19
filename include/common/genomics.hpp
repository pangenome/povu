#ifndef POVU_GENOMIC_TYPES_HPP
#define POVU_GENOMIC_TYPES_HPP

#include <map>
#include <string>
#include <thread>
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

/*setters*/
void append_step(Step s) { this->steps_.emplace_back(s); }
};

typedef Walk StepSeq;

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
const std::vector<Walk> &get_walks() const { return this->walks_; }
std::vector<Walk> &get_walks_mut() { return this->walks_; }

/*setters*/
void set_walks(std::vector<Walk> &&walks) { this->walks_ = walks; }

/* other */
std::string as_str() const {
  return this->fl_.start_.as_str() + this->fl_.end_.as_str();
}
};

 /* a sequence of walks in a RoV --- especially for a given ref */
class It {
  std::vector<Walk> it_;
  pt::idx_t len;

public:
  /* constructor */
  It() : it_(), len(0) {}
  It(Walk &&w) : it_(std::vector<Walk>{w}), len(w.step_count()) {}

  /*getters*/
  pt::idx_t walk_count() const { return this->it_.size(); }
  pt::idx_t step_count() const { return this->len; }

  const std::vector<Walk> &get_walks() const { return this->it_; }
  std::vector<Walk> &get_walks_mut() { return this->it_; }

  const Step &get_step(pt::idx_t idx) const {
    pt::idx_t i = 0;
    for (const Walk &w : this->it_) {
      if (i + w.step_count() > idx) {
        return w.get_steps()[idx - i];
      }
      i += w.step_count();
    }
    throw std::out_of_range("step index out of range");
  }

  const Walk &get_walk_by_step_idx(pt::idx_t idx) const {
    pt::idx_t i = 0;
    for (const Walk &w : this->it_) {
      if (i + w.step_count() > idx) {
        return w;
      }
      i += w.step_count();
    }
    throw std::out_of_range("step index out of range");
  }

  /*setters*/
  void add_walk(Walk &&w) {
    this->it_.emplace_back(w);
    this->len += w.step_count();
  }
};

typedef It Itn;

/*  map of ref_id to the walk of the ref in a RoV */
class RefWalks {
  // use Walk instead of vector<Step>?
  // map of ref_id to the walk of the ref in a RoV
  // a ref can have multiple walks in a RoV
  // the order of walks is not assured
  std::map<pt::id_t, It> ref_walks_;

public:
  /* constructor */
  RefWalks() : ref_walks_() {}

  /*getters*/
  pt::idx_t ref_count() const { return this->ref_walks_.size(); }

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

  /*setters*/
  void add_walk(id_t ref_id, Walk &&w) {
    if (this->ref_walks_.find(ref_id) == this->ref_walks_.end()) {
      this->ref_walks_[ref_id] = It(std::move(w));
    } else {
      this->ref_walks_[ref_id].add_walk(std::move(w));
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
};


/*
  enums
  -----
*/
enum class aln_level_e { step, rov };

} // namespace povu::types::variation

#endif // POVU_GENOMIC_TYPES_HPP
