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
  pt::id_t loop_id_; // the nth time that a ref is going through a flubble RoV
  pt::id_t v_id_;
  pt::idx_t step_idx_;
  pgt::or_e o_;

public:
/* constructor */
Step(pt::id_t v_id, pgt::or_e o)
  : loop_id_(pc::INVALID_ID), v_id_(v_id), step_idx_(pc::INVALID_IDX), o_(o) {}

/*getters*/
pt::idx_t get_step_idx() const { return this->step_idx_; }
pt::id_t get_v_id() const { return this->v_id_; }
pgt::or_e get_o() const { return this->o_; }
pt::id_t get_loop_id() const { return this->loop_id_; }

/*setters*/
void set_step_idx(pt::idx_t step_idx) { this->step_idx_ = step_idx; }
void set_loop_id(pt::id_t loop_id) { this->loop_id_ = loop_id; }
};

class Walk {
std::vector<Step> steps_;

public:
/* constructor */
Walk() : steps_() {}
Walk(pt::id_t id, pgt::or_e o) : steps_(std::vector<Step>{Step{id, o}}) {}

/*getters*/
const std::vector<Step> &get_steps() const { return this->steps_; }

/*setters*/
void append_step(Step s) { this->steps_.emplace_back(s); }
};

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

/*
  typedefs/aliases
  ----------------
*/
// TODO: make these classes?
typedef std::vector<Step> ref_walk;

/*  map of ref_id to the walk of the ref in a RoV */
typedef std::map<id_t, ref_walk> ref_walks;

} // namespace povu::types::variation

#endif // POVU_GENOMIC_TYPES_HPP
