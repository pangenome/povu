#ifndef PV_ALN_HPP
#define PV_ALN_HPP

#include <iostream>
#include <vector>


#include "../common/log.hpp"
#include "../genomics/allele.hpp"

namespace povu::align {
inline constexpr std::string_view MODULE = "povu::align";

namespace pc = povu::constants;
namespace pgt = povu::types::graph;
namespace pga = povu::genomics::allele;

/**
 * used in untangling, the level of alignment
 */

enum class aln_level_e {
  step,
  at // allele traversal
};

struct aln_scores_t {
  pt::idx_t match;
  pt::idx_t mismatch;
  pt::idx_t gap_open;
  pt::idx_t gap_extend;
};

struct aln_result_t {
  pt::idx_t score;
  std::string et; // edit transcript
};

class Matrix {
  std::vector<pt::idx_t> data;
  pt::idx_t row_count_;
  pt::idx_t col_count_;

public:
  Matrix(pt::idx_t row_count, pt::idx_t col_count)
    : data(row_count * col_count, 0), row_count_(row_count),
        col_count_(col_count) {}

  pt::idx_t &operator()(pt::idx_t row, pt::idx_t col) {
    return data[row * this->col_count_ + col];
  }

  pt::idx_t operator()(pt::idx_t row, pt::idx_t col) const {
    return data[row * this->col_count_ + col];
  }

  // Print the matrix to an output stream (default is std::cout)
  void print(std::ostream &os = std::cout) const {
    for (pt::idx_t i = 0; i < this->row_count_; ++i) {
      for (pt::idx_t j = 0; j < this->col_count_; ++j) {
        os << ((*this)(i, j) == pc::INVALID_IDX ? pc::INF : std::to_string((*this)(i, j)))
           << " ";
      }
      os << "\n";
    }
  }
};

std::string align(const pga::itn_t &iw, const pga::itn_t &jw, aln_level_e level);

} // namespace povu::align

#endif // PV_ALN_HPP
