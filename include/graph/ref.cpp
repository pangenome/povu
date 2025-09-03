#include "./ref.hpp"

namespace povu::graph::ref {

bool operator==(const ref_step_t &lhs, const ref_step_t &rhs) {
  return lhs.strand_ == rhs.strand_ && lhs.locus_ == rhs.locus_;
}
bool operator<(const ref_step_t &lhs, const ref_step_t &rhs) {
  if (lhs.strand_ == rhs.strand_) {
    return lhs.locus_ < rhs.locus_;
  }
  return lhs.strand_ < rhs.strand_;
}

} // namespace povu::graph::ref
