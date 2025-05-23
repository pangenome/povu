#ifndef PV_ALN_HPP
#define PV_ALN_HPP

#include <vector>

#include "../common/genomics.hpp"
#include "../common/types.hpp"

namespace povu::align {
namespace pt = povu::types;
namespace pc = povu::constants;
namespace pvt = povu::types::genomics;
namespace pgt = povu::graph_types;

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

std::string align(const pvt::Itn &iw, const pvt::Itn &jw, pvt::aln_level_e level);

} // namespace povu::align

#endif // PV_ALN_HPP
