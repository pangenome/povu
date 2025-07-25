#ifndef POVU_TYPES_GRAPH_HPP
#define POVU_TYPES_GRAPH_HPP

#include <set>

#include "./core.hpp"
#include "./constants.hpp"

namespace povu::types::graph {
namespace pt = povu::types;
namespace pc = povu::constants;

// should this be renamed to clr_e or color_e?
enum class color {
  gray,
  black
};
typedef color color_e;
std::ostream& operator<<(std::ostream& os, const color_e& c);

// Eq class and node id
struct eq_n_id_t {
  std::size_t eq_class;
  std::size_t v_id;
};

struct id_n_cls {
  std::size_t id;
  std::size_t cls;
};


/**
  * l (left) 5' or +
  * r (right) 3' or -
  */
// TODO: replace with struct or class to allow methods like complement
enum class v_end_e {
  l,
  r
};
std::ostream& operator<<(std::ostream& os, const v_end_e& vt);
constexpr v_end_e complement(v_end_e s) {
  return s == v_end_e::l ? v_end_e::r : v_end_e::l;
};

/**
 * l (left) 5' or +
 * r (right) 3' or -
 */
enum class v_type_e {
    l,
    r,
    dummy
};
std::ostream& operator<<(std::ostream& os, const v_type_e& vt);

// TODO: remove
// Merge path_t and biedged PathInfo into one namespace
struct path_t {
  std::string name; // name as pertains the GFA file
  std::size_t id; // numerical id to be associated with handle
  bool is_circular; // is the path circular?
};

struct side_n_id_t {
  v_end_e v_end;
  pt::id_t v_idx; // TODO [B] rename v_idx to v_id

  // -------
  // methods
  // -------
  friend bool operator<(const side_n_id_t& lhs, const side_n_id_t& rhs);
  side_n_id_t complement() const;
};
std::ostream& operator<<(std::ostream& os, const side_n_id_t& x);
typedef side_n_id_t side_n_idx_t; // to use when the id is an index

struct canonical_sese {
  std::size_t start;
  std::size_t end;
  std::set<std::size_t> in_sese; // set of sese ids that contain this sese excluding the start and end
};

enum class or_e {
  forward,
  reverse
};
std::ostream& operator<<(std::ostream& os, const or_e& o);
std::string or_to_str (or_e o);

// TODO: move to povu::types::variation
struct id_or_t {
  pt::id_t v_id; // TODO change type and name to id to pt::id_t
  or_e orientation;

  
  std::string as_str() const {
    return std::format("{}{}", or_to_str(this->orientation) , this->v_id);
  }
};

std::ostream& operator<<(std::ostream& os, const id_or_t& x);
bool operator!=(const id_or_t & lhs, const id_or_t& rhs);
bool operator==(const id_or_t & lhs, const id_or_t& rhs);
bool operator<(const id_or_t& lhs, const id_or_t& rhs);

// a walk is a sequence of vertices also a path
typedef std::vector<id_or_t> walk;


struct flubble {
  id_or_t start_;
  id_or_t end_;

/* constructor(s) */
flubble(id_or_t start, id_or_t end) : start_(start), end_(end) {}

flubble(const std::string& s) {

  // split s based on > and < signs and using s.substr
  // s is in the form of >1<2 or >1>2 or <1<2 or <1>2

  // find the first > or <
  auto first = s.find_first_of("><");
  auto last = s.find_last_of("><");

  // substring based on first and last occurences and store them as size_t

  this->start_.v_id = std::stoull(s.substr(first + 1, last - first - 1));
  this->start_.orientation = s[first] == '>' ? or_e::forward : or_e::reverse;


  this->end_.v_id = std::stoull(s.substr(last + 1, s.size() - last - 1));
  this->end_.orientation = s[last] == '>' ? or_e::forward : or_e::reverse;
}

/* other(s) */
std::string as_str() const {
  std::string s;
  s += this->start_.as_str();
  s += this->end_.as_str();
  return s;
  //return std::format("{}{}", this->start_.as_str(), this->end_.as_str());
}
};

typedef  flubble flubble_t ;

bool operator<(const flubble_t &lhs, const flubble_t &rhs);


class Step {
  // TODO: remove loop_no?
  // pt::id_t loop_no_; // the nth time that a ref is going through a flubble RoV
  pt::id_t v_id_;
  pt::idx_t step_idx_; // also locus
  or_e o_;

public:
  /* constructor */
  Step(pt::id_t v_id, or_e o)
    :v_id_(v_id), step_idx_(pc::INVALID_IDX), o_(o) {}
  Step(pt::id_t v_id, pt::idx_t step_idx, or_e o )
    :v_id_(v_id), step_idx_(step_idx), o_(o) {}

  /*getters*/
  pt::idx_t get_step_idx() const { return this->step_idx_; }
  pt::id_t get_v_id() const { return this->v_id_; }
  or_e get_o() const { return this->o_; }
  //pt::id_t get_loop_no() const { return this->loop_no_; }

  /*setters*/
  void set_step_idx(pt::idx_t step_idx) { this->step_idx_ = step_idx; }
  //void set_loop_no(pt::id_t loop_id) { this->loop_no_ = loop_id; }
};

/* an uninterrupted ordered sequence of steps bound by the start end of a RoV */
class Walk {
protected:
  std::vector<Step> steps_;

public:
  /* constructors */
  Walk() : steps_() {}
  Walk(pt::id_t id, or_e o) : steps_(std::vector<Step>{Step{id, o}}) {}
  Walk(Step s) : steps_(std::vector<Step>{s}) {}

  /*getters*/
  pt::idx_t step_count() const { return this->steps_.size(); }
  const std::vector<Step> &get_steps() const { return this->steps_; }
  std::vector<Step> &get_steps_mut() { return this->steps_; }
  const Step &get_step(pt::idx_t idx) const { return this->steps_[idx]; }
  std::string as_str() const {
    std::string s;
    for (const Step &step : this->steps_) {

      s += std::format("{}{}", step.get_o() == or_e::forward ? ">" : "<",
                       step.get_v_id());
    }
    return s;
  }

  /*setters*/
  void append_step(Step s) { this->steps_.emplace_back(s); }
};

typedef Walk StepSeq;


} // namespace povu::graph_types

#endif
