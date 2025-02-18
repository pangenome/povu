#ifndef PV_ALN_HPP
#define PV_ALN_HPP

#include <vector>

#include "../common/genomics.hpp"
#include "../common/types.hpp"

namespace povu::align {
namespace pt = povu::types;
namespace pc = povu::constants;
namespace pvt = povu::types::genomics;

struct aln_scores_t {
  pt::idx_t match;
  pt::idx_t mismatch;
  pt::idx_t gap_open;
  pt::idx_t gap_extend;
};

struct aln_result_t {
  pt::idx_t score;
  std::string alignment;
};

void align_steps(const pvt::ref_walk &iw, const pvt::ref_walk &jw);
void align_rovs(const pvt::ref_walk &iw, const pvt::ref_walk &jw);

} // namespace povu::align

#endif // PV_ALN_HPP
