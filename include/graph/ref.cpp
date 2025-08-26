#include "./ref.hpp"

namespace povu::graph::ref {
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
} // namespace povu::graph::ref
