#ifndef POVU_GENOMIC_TYPES_HPP
#define POVU_GENOMIC_TYPES_HPP

#include <map>
#include <memory>
#include <string>
#include <thread>
#include <utility>
#include <vector>

//#include "../../graph/bidirected.hpp"
#include "./pvst.hpp"
#include "./types.hpp"
#include "./graph.hpp"

namespace povu::types::genomics {
namespace pgt = povu::types::graph;
namespace pt = povu::types;
namespace pc = povu::constants;
namespace pvst = povu::types::pvst;



/* region of variation */
class RoV {
std::vector<pgt::Walk> walks_;
const pvst::VertexBase *pvst_vtx;

public:

// --------------
// constructor(s)
// --------------

RoV(const pvst::VertexBase *v) : walks_(), pvst_vtx(v) {}

// --------------
// getter(s)
// --------------

pt::idx_t walk_count() const { return this->walks_.size(); }

// const pgt::id_or_t &get_entry() const { return this->fl_.start_; }
// const pgt::id_or_t &get_exit() const { return this->fl_.end_; }
// pgt::flubble_t get_flb() const { return this->fl_; }
const std::vector<pgt::Walk> &get_walks() const { return this->walks_; }
std::vector<pgt::Walk> &get_walks_mut() { return this->walks_; }

// --------------
// setter(s)
// --------------

void set_walks(std::vector<pgt::Walk> &&walks) { this->walks_ = walks; }

// --------------
// other(s)
// --------------

std::string as_str() const {
  return this->pvst_vtx->as_str();
}
};

/* variation related types */

/*allele traversal---a walk taken by a reference*/
// when the AT is a deletion, the walk is empty
//typedef Walk AT;
class AT : public pgt::Walk {
bool is_del_;
public:
  /* constructors */
  AT() : pgt::Walk(){}
  AT(pgt::Step s) : pgt::Walk(s) {}

  /*getters*/
  bool is_del() const { return this->is_del_; }

  /* setters */
  void set_is_del(bool is_del) { this->is_del_ = is_del; }
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
  //It(AT &&w) : it_(std::vector<AT>{w}), len(w.step_count()) {}

  /* getters */
  pt::idx_t at_count() const { return this->it_.size(); }
  pt::idx_t step_count() const { return this->len; }

  const std::vector<AT> &get_ats() const { return this->it_; }
  std::vector<AT> &get_ats_mut() { return this->it_; }

  const AT &get_at(pt::idx_t at_idx) const { return this->it_[at_idx]; }
  AT &get_at_mut(pt::idx_t at_idx) { return this->it_[at_idx]; }

  const pgt::Step &get_step(pt::idx_t idx) const {
    pt::idx_t i = 0;
    for (const AT &w : this->it_) {
      if (i + w.step_count() > idx) {
        return w.get_steps()[idx - i];
      }
      i += w.step_count();
    }
    throw std::out_of_range("step index out of range");
  }

  // const AT &get_at_by_step_idx(pt::idx_t idx) const {
  //   pt::idx_t i = 0;
  //   for (const AT &w : this->it_) {
  //     if (i + w.step_count() > idx) {
  //       return w;
  //     }
  //     i += w.step_count();
  //   }
  //   throw std::out_of_range("step index out of range");
  // }

  /* setters */
  void append_at(AT &&w) {
    this->it_.emplace_back(w);
    this->len += w.step_count();
  }

  void append_step(pt::idx_t at_idx, pgt::Step s) {
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

  //void dec_step_count(pt::idx_t dec) { this->len -= dec; }
};
typedef It Itn;

/*  map of ref_id to the walk of the ref in a RoV */
class RefWalks {

  // use Walk instead of vector<Step>?
  // map of ref_id to the walk of the ref in a RoV
  // a ref can have multiple walks in a RoV
  // the order of walks is not assured
  std::map<pt::id_t, Itn> ref_itns_;

  //alignment between two refs
  std::map<pt::up_t<pt::id_t>, std::string> aln;

  pgt::flubble_t fl_;

  bool is_tangled_;

  /* private methods */
  // returns an unordered pair
  // up_t to_up (pt::id_t a, pt::id_t b) const {
  //   return {std::min(a, b), std::max(a, b)};
  // };

public:
  /* constructor */
  RefWalks(pgt::flubble_t fl) : ref_itns_(), fl_(fl), is_tangled_(false) {}

  /*getters*/
  pt::idx_t ref_count() const { return this->ref_itns_.size(); }
  const pgt::flubble_t &get_flb() const { return this->fl_; }

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

  const std::map<pt::id_t, Itn> &get_ref_itns() const {
    return this->ref_itns_;
  }

  std::map<pt::id_t, Itn> &get_ref_itns_mut() {
    return this->ref_itns_;
  }

  // bool has_aln(pt::id_t ref_id1, pt::id_t ref_id2) const {
  //   return this->aln.find(to_up(ref_id1, ref_id2)) != this->aln.end();
  // }

  const std::string &get_aln(pt::id_t ref_id1, pt::id_t ref_id2) const {
    return this->aln.at(pt::up_t<pt::id_t>{ref_id1, ref_id2});
  }

  const std::map<pt::up_t<pt::id_t>, std::string> &get_alns() const {
    return this->aln;
  }

  bool is_tangled() const { return this->is_tangled_; }

  /*setters*/
  // returns true to mean that tangling exists. a walk traverses an RoV more than once
  void add_itn(id_t ref_id, Itn &&itn) {
    if (this->ref_itns_.find(ref_id) == this->ref_itns_.end()) {
      this->ref_itns_[ref_id] = std::move(itn);
    } else {
      Itn &itn_ = this->ref_itns_[ref_id];
      for (auto &w : itn.get_ats_mut()) {
        itn_.append_at(std::move(w));
      }

      //is_tangled_ = true;
    }
  }

  void replace_itn(pt::id_t ref_id, Itn &&itn) {
    this->ref_itns_[ref_id] = std::move(itn);
  }

  void sort_by_step_idx() {

    auto compare_si = [](const AT &a, const AT &b) {
      return a.get_steps().front().get_step_idx() < b.get_steps().front().get_step_idx();
    };

    // for each ref, sort the walks by step_idx of the first step
    for (auto &[_, itn] : this->ref_itns_) {
      std::sort(itn.get_ats_mut().begin(), itn.get_ats_mut().end(), compare_si);
    }
  }

    // set the loop number for each step
    // steps within the same walk share the same loop number
    // for (auto &[_, itn] : this->ref_walks_) {
    //   //pt::id_t loop_no = 0;
    //   for (Walk &w : itn.get_walks_mut()) {
    //     for (Step &s : w.get_steps_mut()) {
    //       // s.set_loop_no(loop_no);
    //     }
    //     //loop_no++;
    //   }

  void add_aln(pt::id_t ref_id1, pt::id_t ref_id2 , std::string &&aln) {
    this->aln[pt::up_t<pt::id_t>{ref_id1, ref_id2}] = aln;
  }

  void set_tangled(bool is_tangled) { this->is_tangled_ = is_tangled; }
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
  VcfRec(pt::id_t ref, pt::idx_t pos, std::string id, AT ref_at, std::vector<AT> alt_ats)
      : ref(ref), pos(pos), id(id), ref_at(ref_at), alt_ats(alt_ats) {}

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
