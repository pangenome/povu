#include <iostream>

#include "./types.hpp"

namespace povu::graph_types {

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

/*
 * Orientation
 * ------------
 */
// >> and << might be better than + and -
std::ostream& operator<<(std::ostream& os, const orientation_t& o) {
  switch (o) {
  case orientation_t::forward:
  os << ">";
  break;
  case orientation_t::reverse:
  os << "<";
  break;
  }

  return os;
}
std::string or_to_str (orientation_t o) {
  return o == orientation_t::forward ? ">" : "<";
};

/*
 * id and orientation
 * -----------------
 */
std::ostream& operator<<(std::ostream& os, const id_n_orientation_t& x) {
  os << x.orientation << x.v_idx;
  return os;
}

bool operator<(const id_n_orientation_t& lhs, const id_n_orientation_t& rhs) {
  if (lhs.v_idx < rhs.v_idx) {
    return true;
  }
  else if (lhs.v_idx == rhs.v_idx) {
    return lhs.orientation < rhs.orientation;
  }
  else {
    return false;
  }
}

bool operator!=(const id_n_orientation_t & lhs, const id_n_orientation_t& rhs) {
    return lhs.v_idx != rhs.v_idx || lhs.orientation != rhs.orientation;
}

bool operator==(const id_n_orientation_t & lhs, const id_n_orientation_t& rhs) {
    return lhs.v_idx == rhs.v_idx && lhs.orientation == rhs.orientation;
}
} // namespace povu::graph_types
