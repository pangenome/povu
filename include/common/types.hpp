#ifndef POVU_TYPES_HPP
#define POVU_TYPES_HPP

#include <chrono>
#include <cstddef>
#include <string>
#include <utility>
#include <format>
#include <iostream>
#include <set>
#include <sys/types.h>
#include <utility>
#include <vector>

namespace povu::types {

typedef std::chrono::high_resolution_clock Time; // C++ timer

typedef u_int32_t id_t;
typedef u_int32_t idx_t;

struct Stride {
  std::size_t start;
  std::size_t length;
};
typedef Stride span;


typedef std::pair<std::size_t, std::size_t> size_t_pair;

/**
 * ordered pair similar to std::pair but with same type on both sides for less
 * typing
 */
template <typename T> struct Pair {
  T first;
  T second;

  friend bool operator<(const Pair &lhs, const Pair &rhs) {
    return std::tie(lhs.first, lhs.second) < std::tie(rhs.first, rhs.second);
  }

  // spaceship operator
  friend constexpr auto operator<=>(Pair, Pair) = default;
};

template <typename T> struct unordered_pair {
  T l;
  T r;

  unordered_pair(T l, T r) : l(std::min(l, r)), r(std::max(l, r)) {}

  // spaceship operator
  friend constexpr auto operator<=>(const unordered_pair &,
                                    const unordered_pair &) = default;
};

} // namespace povu::types

namespace povu::constants {
//
const std::size_t DUMMY_VERTEX_COUNT { 2 };

// colors
const std::string gray{"gray"};
const std::string black{"black"};
const std::string red{"red"};

// numeric
const std::size_t SIZE_T_MIN = std::numeric_limits<size_t>::min();
const std::size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();
const int UNDEFINED_INT = std::numeric_limits<int>::min();
const std::size_t UNDEFINED_SIZE_T = std::numeric_limits<size_t>::max();
const std::size_t UNDEFINED_IDX = std::numeric_limits<povu::types::idx_t>::max();
const std::size_t UNDEFINED_ID = std::numeric_limits<povu::types::id_t>::max();
const std::size_t DUMMY_VTX_ID = UNDEFINED_ID;
const std::size_t INVALID_ID = UNDEFINED_SIZE_T;
const std::size_t INVALID_IDX = UNDEFINED_SIZE_T;

// strings
const std::string EMPTY_SET = "\u2205";
const std::string UNDEFINED_VALUE = "\u2205";
const std::string WAVY_ARROW = "\u2933";

// genomics constants
const std::string UNDEFINED_PATH_LABEL { "undefined" };
const std::size_t UNDEFINED_PATH_ID { INVALID_ID };
const std::size_t UNDEFINED_PATH_POS { INVALID_ID };

// VCF
const char COL_SEP = '\t'; // column separator
const char NO_VALUE = '.'; // null character
} // namespace povu::constants



namespace povu::graph_types {

/*
 * black edge is default
 * TODO: pick a better default
 * gray edge is a bi-edge
 */
enum class color { gray, black };
typedef  color color_e;
std::ostream& operator<<(std::ostream& os, const color& c);

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
enum class VertexEnd {
  l,
  r
};
typedef VertexEnd v_end_t;
typedef VertexEnd v_end;
typedef VertexEnd VtxEnd;

std::ostream& operator<<(std::ostream& os, const VertexEnd& vt);
VertexEnd complement(VertexEnd s);


/**
 * l (left) 5' or +
 * r (right) 3' or -
 */
enum class VertexType {
    l,
    r,
    dummy
};
typedef VertexType v_type; // deprected use v_type_e
typedef VertexType v_type_e;
std::ostream& operator<<(std::ostream& os, const VertexType& vt);

// Merge path_t and biedged PathInfo into one namespace
struct path_t {
  std::string name; // name as pertains the GFA file
  std::size_t id; // numerical id to be associated with handle
  bool is_circular; // is the path circular?
};

struct side_n_id_t {
  VertexEnd v_end;
  std::size_t v_idx;

  // -------
  // methods
  // -------
  friend bool operator<(const side_n_id_t& lhs, const side_n_id_t& rhs);
  side_n_id_t complement() const;
};
std::ostream& operator<<(std::ostream& os, const side_n_id_t& x);

struct canonical_sese {
  std::size_t start;
  std::size_t end;
  std::set<std::size_t> in_sese; // set of sese ids that contain this sese excluding the start and end
};

enum class orientation_t {
  forward,
  reverse
};
typedef orientation_t or_t;
std::ostream& operator<<(std::ostream& os, const orientation_t& o);
std::string or_to_str (orientation_t o);

// TODO rename to id_or ?
struct id_n_orientation_t {
  std::size_t v_idx; // rename to id
  orientation_t orientation;

  std::string as_str() const {
    return std::format("{}{}", or_to_str(this->orientation) , this->v_idx);
  }
};
typedef id_n_orientation_t id_or;

std::ostream& operator<<(std::ostream& os, const id_n_orientation_t& x);
bool operator!=(const id_n_orientation_t & lhs, const id_n_orientation_t& rhs);
bool operator==(const id_n_orientation_t & lhs, const id_n_orientation_t& rhs);
bool operator<(const id_n_orientation_t& lhs, const id_n_orientation_t& rhs);

typedef std::vector<graph_types::id_n_orientation_t> walk; // a walk is a sequence of vertices also a path

struct flubble {
  id_n_orientation_t start_;
  id_n_orientation_t end_;

// constructors

flubble(id_or start, id_or end) : start_(start), end_(end) {}

flubble(const std::string& s) {

  // split s based on > and < signs and using s.substr
  // s is in the form of >1<2 or >1>2 or <1<2 or <1>2

  // find the first > or <
  auto first = s.find_first_of("><");
  auto last = s.find_last_of("><");

  // substring based on first and last occurences and store them as size_t

  this->start_.v_idx = std::stoull(s.substr(first + 1, last - first - 1));
  this->start_.orientation = s[first] == '>' ? orientation_t::forward : orientation_t::reverse;


  this->end_.v_idx = std::stoull(s.substr(last + 1, s.size() - last - 1));
  this->end_.orientation = s[last] == '>' ? orientation_t::forward : orientation_t::reverse;

}
};

} // namespace povu::graph_types
#endif
