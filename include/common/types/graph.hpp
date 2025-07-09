#ifndef POVU_TYPES_GRAPH_HPP
#define POVU_TYPES_GRAPH_HPP

#include <set>

#include "./core.hpp"

namespace povu::types::graph {
namespace pt = povu::types;

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
} // namespace povu::graph_types

#endif
