#include "./align.hpp"


namespace povu::align {

std::vector<std::vector<pt::idx_t>> matrix;

class Matrix {
  std::vector<pt::idx_t> data;
  // pt::idx_t row_count_;
  pt::idx_t col_count_;

public:
  Matrix(pt::idx_t row_count, pt::idx_t col_count)
    : data(row_count * col_count, 0), col_count_() {}

  pt::idx_t &operator()(pt::idx_t row, pt::idx_t col) {
    return data[row * col_count_ + col];
  }
};

inline pt::idx_t min(pt::idx_t a, pt::idx_t b) {
  return std::min(a, b);
}
inline pt::idx_t min(pt::idx_t a, pt::idx_t b, pt::idx_t c) {
  return std::min(a, std::min(b, c));
}


// smith waterman alignment, affine gap penalty
aln_result_t sw_align(pt::idx_t col_count, pt::idx_t row_count,
                      const std::string &str1, const std::string &str2,
                      bool (*eq)(pt::idx_t, pt::idx_t, const std::string &str1,
                                 const std::string &str2)) {

  aln_scores_t scores = {
    .match = 0,
    .mismatch = 1,
    .gap_open = 2,
    .gap_extend = 1
  };

  // init M , I , D matrices
  Matrix M(row_count, col_count);
  Matrix I(row_count, col_count);
  Matrix D(row_count, col_count);

  // init I top row
  for (pt::idx_t i = 0; i < col_count; ++i) {
    I(0, i) = pc::INVALID_IDX;
  }

  // init D left column
  for (pt::idx_t i = 0; i < row_count; ++i) {
    D(i, 0) = pc::INVALID_IDX;
  }

  // init M
  M(0, 0) = 0;

  // fill M, I, D matrices
  for (pt::idx_t i = 1; i < row_count; ++i) { // rows
    for (pt::idx_t j = 1; j < col_count; ++j) { // cols
      // fill M
      M(i, j) = min(M(i - 1, j - 1) + (eq(i,j, str1, str2) ? scores.match : scores.mismatch),
                    I(i - 1, j),
                    D(i, j - 1));

      // fill I
      I(i, j) = min(
        M(i - 1, j) + scores.gap_open + scores.gap_extend,
        I(i - 1, j) + scores.gap_extend
      );

      // fill D
      D(i, j) = min(
        M(i, j - 1) + scores.gap_open + scores.gap_extend,
        D(i, j - 1) + scores.gap_extend
      );

    }
  }

  pt::idx_t aln_score = M(row_count - 1, col_count - 1);
  std::cout << "aln score: " << aln_score << std::endl;

  // backtrace
  std::string aln;
  pt::idx_t i = row_count - 1;
  pt::idx_t j = col_count - 1;
  while (i > 0 && j > 0) {
    pt::idx_t score = M(i, j);
    pt::idx_t score_diag = M(i - 1, j - 1);
    pt::idx_t score_up = I(i, j);
    pt::idx_t score_left = D(i, j);

    if (score == score_diag + (eq(i,j, str1, str2) ? scores.match : scores.mismatch)) {
      aln.push_back('M');
      --i;
      --j;
    } else if (score == score_up) {
      aln.push_back('I');
      --i;
    } else if (score == score_left) {
      aln.push_back('D');
      --j;
    } else {
      std::cerr << "error" << std::endl;
      break;
    }
  }

  std::reverse(aln.begin(), aln.end());

  return {aln_score, aln};
}

bool eq_step(pt::idx_t x, pt::idx_t y) {}

bool eq_rov(pt::idx_t x, pt::idx_t y) {}

void align_steps(const pvt::ref_walk &iw, const pvt::ref_walk &jw) {

  // two steps are equal if they
  // share the same vertex id, the same loop id, and the same orientation
  // false otherwise
  pt::idx_t i_step_count = iw.size();
  pt::idx_t j_step_count = jw.size();

  sw_align(i_step_count, j_step_count, eq_step);

}


void align_rovs(const pvt::ref_walk &iw, const pvt::ref_walk &jw) {


}

bool eq_chars(pt::idx_t idx1,
              pt::idx_t idx2,
              const std::string &str1,
              const std::string &str2) {
  return str1[idx1] == str2[idx2];
}
void align_strings(const std::string &str1, const std::string &str2) {
  sw_align(str1.length(), str2.length(), str1, str2, eq_chars);
}

} // namespace povu::align
