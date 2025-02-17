#include "./types.hpp"
#include <string>
#include <vector>
#include <thread>
#include <map>

namespace povu::types::variation {
namespace pgt = povu::graph_types;
namespace pt = povu::types;

struct ref_step_t {
  pt::id_t loop_id_; // the nth time that a ref is going through a flubble RoV
  pt::id_t v_id_;
  pt::idx_t step_idx;
  pgt::or_e o_;

/* constructor */
ref_step_t(pt::id_t loop_id, pt::id_t v_id, pt::idx_t step_idx, pgt::or_e o)
  : loop_id_(loop_id), v_id_(v_id), step_idx(step_idx), o_(o) {}

/* getter(s) */
pt::id_t get_v_id() const { return this->v_id_; }

/* setter(s) */
void set_loop_id(pt::id_t loop_id) { this->loop_id_ = loop_id; }
};

// vector of ref_steps
typedef std::vector<ref_step_t> ref_walk;
// map of ref_id to vector of ref_steps
typedef std::map<id_t, ref_walk> ref_walks;

typedef pgt::id_or_t step;

class Walk {
std::vector<step> steps_;

public:
/* constructor */
Walk() : steps_() {}
Walk(pt::id_t id, pgt::or_e o) : steps_(std::vector<step>{step{id, o}}) {}

/*getters*/
const std::vector<step> &get_steps() const { return this->steps_; }

/*setters*/
void append_step(step s) { this->steps_.emplace_back(s); }
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

} // namespace povu::types::variation
