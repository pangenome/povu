#ifndef POVU_TYPES_HPP
#define POVU_TYPES_HPP

#include <chrono>
#include <cstddef>
#include <cstdint>
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
typedef int8_t status_t; // return status of a fn

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

/**
 * unordered pair with same type on both sides
 */
template <typename T> struct unordered_pair {
  T l;
  T r;

  unordered_pair(T l, T r) : l(std::min(l, r)), r(std::max(l, r)) {}

  // spaceship operator
  friend constexpr auto operator<=>(const unordered_pair &,
                                    const unordered_pair &) = default;
};

template <typename T>
using up_t = unordered_pair<T>;
// up = unordered pair
//typedef std::pair<id_t, id_t> up_t;

} // namespace povu::types


namespace povu::graph_types {
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

namespace povu::constants {
namespace pt = povu::types;

// colors
const std::string GRAY{"gray"};
const std::string BLACK{"black"};
const std::string RED{"red"};
const std::string BLUE{"blue"};

// numeric
const std::size_t SIZE_T_MIN = std::numeric_limits<size_t>::min();
const std::size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();
const int UNDEFINED_INT = std::numeric_limits<int>::min();
const std::size_t UNDEFINED_SIZE_T = SIZE_T_MAX;
const pt::idx_t MAX_ID = std::numeric_limits<povu::types::idx_t>::max();
const pt::idx_t MAX_IDX = std::numeric_limits<povu::types::idx_t>::max();

const pt::idx_t UNDEFINED_IDX = MAX_IDX; // TODO: replace with INVALID
const pt::id_t UNDEFINED_ID = MAX_ID;    // TODO: replace with INVALID
const pt::id_t DUMMY_VTX_ID = UNDEFINED_ID;
const pt::id_t INVALID_ID = MAX_ID;
const pt::idx_t INVALID_IDX = MAX_IDX;
const pt::idx_t INVALID_CLS = MAX_IDX; // equivalence class

// strings
const std::string EMPTY_SET = "\u2205";
const std::string UNDEFINED_VALUE = EMPTY_SET;
const std::string WAVY_ARROW = "\u2933";
const std::string INF = "\u221E"; //infinity

// genomics constants
const std::string UNDEFINED_PATH_LABEL{"undefined"};
const std::size_t UNDEFINED_PATH_ID{INVALID_ID};
const std::size_t UNDEFINED_PATH_POS{INVALID_ID};

const char PVST_HEADER_SYMBOL = 'H';
const char PVST_FLUBBLE_SYMBOL = 'F';

// VCF
const char COL_SEP = '\t'; // column separator
const char NO_VALUE = '.'; // null character
} // namespace povu::constants



#endif
