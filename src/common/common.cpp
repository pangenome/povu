#include <cstddef>
#include <iostream>
#include <utility>

#include "./common.hpp"

namespace common_fns {

std::size_t to_bidirected_idx(std::size_t x, bool has_dummy) {
  // we added 1 because we added a dummy start node before bi-edging
  if (has_dummy) { --x; }
  return x % 2 == 0 ? ((x + 2) / 2) - 1 : ((x + 1) / 2) - 1;
}


std::pair<std::size_t, std::size_t> frm_bidirected_idx(std::size_t x, bool has_dummy) {
  auto foo = [&](std::size_t i) { return (2*(i+1)) + (has_dummy ?  1 : 0 ); };
  return { foo(x) - 2, foo(x) - 1 };
}

} // namespace common_fns

namespace graph_types {

std::ostream& operator<<(std::ostream& os, const VertexType& vt) {
  switch (vt) {
  case VertexType::l:
    os << "+";
    break;
  case VertexType::r:
    os << "-";
    break;
  default:
    os << "*";
    break;
  }

  return os;
}

/*
 * VertexEnd
 * ---------
 */
std::ostream& operator<<(std::ostream& os, const VertexEnd& ve) {
  switch (ve) {
  case VertexEnd::l:
  os << "+";
  break;
  case VertexEnd::r:
  os << "-";
  break;
  }

  return os;
}

VertexEnd complement(VertexEnd s) {
  return s == VertexEnd::l ? VertexEnd::r : VertexEnd::l;
};


std::ostream& operator<<(std::ostream& os, const color& c) {
  switch (c) {
  case color::gray:
    os << "gray";
    break;
  case color::black:
    os << "black";
    break;
  default:
    os << "unknown";
    break;
  }
  return os;
}


/*
 * Side and SideID
 * ----
 */
bool operator<(const side_n_id_t& lhs, const side_n_id_t& rhs) {
  if (lhs.v_idx < rhs.v_idx) {
    return true;
  }
  else if (lhs.v_idx == rhs.v_idx) {
    return lhs.v_end < rhs.v_end;
  }
  else {
    return false;
  }
}

std::ostream& operator<<(std::ostream& os, const side_n_id_t& x) {
  os << x.v_idx << x.v_end;
  return os;
}

side_n_id_t side_n_id_t::complement() const {
  return {graph_types::complement(this->v_end), this->v_idx};
}
} // namespace graph_types
