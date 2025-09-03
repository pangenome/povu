#ifndef POVU_TYPES_GRAPH_HPP
#define POVU_TYPES_GRAPH_HPP

#include <iostream>
#include <string>
#include <string_view>
#include <vector>

#include "../common/core.hpp"
#include "../common/compat.hpp"


namespace povu::types::graph {

// should this be renamed to clr_e or color_e?
enum class color_e {
  gray,
  black
};
std::string_view to_str(color_e ve);
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
std::string_view to_str(v_end_e ve);
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
std::string_view to_str(v_type_e vt);

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


enum class or_e { forward, reverse };
std::string_view to_str (or_e o);
std::ostream& operator<<(std::ostream& os, const or_e& o);


struct id_or_t {
  pt::id_t v_id;
  or_e orientation;

  std::string as_str() const {
    return pv_cmp::format("{}{}", to_str(this->orientation) , this->v_id);
  }
};

std::ostream& operator<<(std::ostream& os, const id_or_t& x);
bool operator!=(const id_or_t & lhs, const id_or_t& rhs);
bool operator==(const id_or_t & lhs, const id_or_t& rhs);
bool operator<(const id_or_t& lhs, const id_or_t& rhs);

typedef id_or_t step_t; // TODO: rename to graph_step_t
typedef std::vector<id_or_t> walk_t;



} // namespace povu::graph_types

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace ptg = povu::types::graph;

#endif
