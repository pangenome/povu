#include "./graph.hpp"

namespace povu::types::graph {

std::ostream& operator<<(std::ostream& os, const v_type_e& vt) {
  switch (vt) {
  case v_type_e::l:
    os << "+";
    break;
  case v_type_e::r:
    os << "-";
    break;
  default:
    os << "*";
    break;
  }

  return os;
}

/*
 * v_end_e
 * ---------
 */
std::ostream& operator<<(std::ostream& os, const v_end_e& ve) {
  switch (ve) {
  case v_end_e::l:
  os << "+";
  break;
  case v_end_e::r:
  os << "-";
  break;
  }

  return os;
}

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
  return {types::graph::complement(this->v_end), this->v_idx};
}

/*
 * Orientation
 * ------------
 */
// >> and << might be better than + and -
std::ostream& operator<<(std::ostream& os, const or_e& o) {
  switch (o) {
  case or_e::forward:
  os << ">";
  break;
  case or_e::reverse:
  os << "<";
  break;
  }

  return os;
}
std::string or_to_str (or_e o) {
  return o == or_e::forward ? ">" : "<";
};

/*
 * id and orientation
 * -----------------
 */
std::ostream& operator<<(std::ostream& os, const id_or_t& x) {
  os << x.orientation << x.v_id;
  return os;
}

bool operator<(const id_or_t& lhs, const id_or_t& rhs) {
  if (lhs.v_id < rhs.v_id) {
    return true;
  }
  else if (lhs.v_id == rhs.v_id) {
    return lhs.orientation < rhs.orientation;
  }
  else {
    return false;
  }
}

bool operator!=(const id_or_t & lhs, const id_or_t& rhs) {
    return lhs.v_id != rhs.v_id || lhs.orientation != rhs.orientation;
}

bool operator==(const id_or_t & lhs, const id_or_t& rhs) {
    return lhs.v_id == rhs.v_id && lhs.orientation == rhs.orientation;
}



std::optional<pan_sn> label_to_pan_sn(const std::string &label, char delim) {

  pan_sn res;

  // split on delim
  std::vector<std::string> parts;
  {
    std::string token;
    std::istringstream in{label};
    while (std::getline(in, token, delim)) {
      parts.push_back(std::move(token));
    }
  }

  // we expect exactly 3 parts : sample, haplotype, contig
  if (parts.size() == 3) {
    auto &s = parts[0];
    auto &h = parts[1]; // we expect the haplotype to be a number
    auto &c = parts[2];

    // haplotype must be all digits
    if (pu::is_numeric_string(h)) {
      res.sample_name_ = std::move(s);
      res.haplotype_id_ = static_cast<pt::id_t>(std::stoull(h));
      res.contig_name_ = std::move(c);
    }
    // if not all digits, return empty optional
    else {
      return std::nullopt;
    }
  } else {
    // if not 3 parts, return empty optional
    return std::nullopt;
  }

  return std::optional<pan_sn>{res};
}

} // namespace povu::graph_types
