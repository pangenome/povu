#include "./align.hpp"

namespace povu::align {
// match result
struct match_res_t {
  pt::idx_t a_inc;
  pt::idx_t b_inc;
  bool is_match;
};
inline pt::idx_t min(pt::idx_t a, pt::idx_t b) { return std::min(a, b); }
inline pt::idx_t min(pt::idx_t a, pt::idx_t b, pt::idx_t c) {
  return std::min(a, std::min(b, c));
}

inline match_res_t eq_step(const pvt::Itn &a, pt::idx_t a_idx,
                           const pvt::Itn &b, pt::idx_t b_idx) {

  if (a_idx >= a.step_count() || b_idx >= b.step_count()) {
    std::cerr << "Index out of bounds\n";
    exit(1);
  }

  auto is_match = [](const pgt::Step &a, const pgt::Step &b) {
    return a.get_v_id() == b.get_v_id() && a.get_o() == b.get_o();
  };


  const pgt::Step &a_step = a.get_step(a_idx);
  const pgt::Step b_step = b.get_step(b_idx);


  if (is_match(a_step, b_step)) {
    return {1, 1, true};
  }

  return {1, 1, false};
}

/*for RoV inc by the length of the walk*/
inline match_res_t eq_at(const pvt::Itn &a, pt::idx_t a_idx,
                          const pvt::Itn &b, pt::idx_t b_idx) {

  // if any of the steps is not a match in the ROV then it is not a match
  const pvt::AT &a_at = a.get_at(a_idx);
  const pvt::AT &b_at = b.get_at(b_idx);


  pt::idx_t a_jmp = a_at.step_count();
  pt::idx_t b_jmp = b_at.step_count();

  if (a_jmp != b_jmp) {
    return {1, 1, false};
  }

  // TODO: also compare loop no
  auto is_match = [](const pgt::Step &a, const pgt::Step &b) {
    return a.get_v_id() == b.get_v_id() && a.get_o() == b.get_o();
  };

  // check for the order as well

  for (pt::idx_t i {}; i < a_jmp; i++) {
    if (!is_match(a_at.get_step(i), b_at.get_step(i))) {
      return {1, 1, false};
    }
  }

  return {1, 1, true};
}

aln_result_t global_align(const pvt::Itn &str1, pt::idx_t str1_len,
                          const pvt::Itn &str2, pt::idx_t str2_len,
                          const aln_scores_t &scores,
                          match_res_t (*eq)(const pvt::Itn &, pt::idx_t, const pvt::Itn &, pt::idx_t)) {
  // Define the scoring parameters
  const pt::idx_t a = scores.match;
  const pt::idx_t x = scores.mismatch;
  const pt::idx_t o = scores.gap_open;
  const pt::idx_t e = scores.gap_extend;

  pt::idx_t row_count = str1_len + 1;
  pt::idx_t col_count = str2_len + 1;

  // Initialize M, I, D matrices.
  Matrix M(row_count, col_count);
  Matrix I(row_count, col_count);
  Matrix D(row_count, col_count);

  // Initialize the top row of I.
  for (pt::idx_t j = 1; j < col_count; ++j) {
    D(0, j) = pc::INVALID_IDX;
    I(0, j) = (j * e) + o;
  }

  // Initialize the left column of D.
  for (pt::idx_t i = 1; i < row_count; ++i) {
    I(i, 0) = pc::INVALID_IDX;
    D(i, 0) = (i * e) + o;
  }

  // Initialize M(0,0)
  M(0, 0) = 0;
  D(0,0) = pc::INVALID_IDX;
  I(0,0) = pc::INVALID_IDX;
  // Initialize first column: gap penalties for aligning str1 with an empty string.
  for (pt::idx_t i = 1; i < row_count; ++i) {
    M(i, 0) = pc::INVALID_IDX;
  }
  // Initialize first row.
  for (pt::idx_t j = 1; j < col_count; ++j) {
    M(0, j) = pc::INVALID_IDX;
  }


  // Fill M, I, D matrices.
  // TODO: is it possible to depend on the i and j inc values?
  for (pt::idx_t i = 1; i < row_count; i++) {  // rows
    for (pt::idx_t j = 1; j < col_count; j++) {  // cols

      // Fill I (gap in str2, insertion in str1)
      I(i, j) = min(
        M(i - 1, j) + o + e,
        I(i - 1, j) + e
      );

      // Fill D (gap in str1, deletion in str1)
      D(i, j) = min(
        M(i, j - 1) + o + e,
        D(i, j - 1) + e
      );

      // Fill M (match/mismatch from diagonal).
      // Use i-1 and j-1 for the characters from the strings.
      auto [_, __, is_match] = eq(str1, i - 1, str2, j - 1);

      //std::cerr << "i: " << i << " j: " << j << " is_match: " << is_match << "\n";

      M(i, j) = min(
        M(i - 1, j - 1) + (is_match ? a : x),
        I(i, j),
        D(i, j)
      );
    }
  }

  pt::idx_t aln_score = M(row_count-1, col_count - 1);


  //std::cerr << "Score: " << aln_score << "\n";

  // Define an enum to track which matrix we’re in.


  pt::idx_t i = str1_len;
  pt::idx_t j = str2_len;
  std::string et;
  et.reserve(std::max(str1_len, str2_len));

  if (false) { // Debug output: print matrices.
    std::cerr << "[" << row_count << "," << col_count << "]\n";
    std::cerr << "M Matrix:\n";
    M.print(std::cerr);
    std::cerr << "I Matrix:\n";
    I.print(std::cerr);
    std::cerr << "D Matrix:\n";
    D.print(std::cerr);
    std::cerr << "\n";
  }

  /* Trace back to reconstruct the alignment. */
  while (i > 0 && j > 0) {
    auto [dec_i, dec_j, is_match] = eq(str1, i - 1, str2, j - 1);
    // std::cerr << "i: " << i << " j: " << j << " is_match: " << is_match << "\n";

    if (M(i,j) == min(M(i-1,j), M(i, j-1), (M(i-1, j-1) + (is_match ? a : x)))) {
      //std::cerr << "ps M\n";
      et.push_back((is_match ? 'M' : 'X'));
      i -= 1;
      j -= 1;
    }
    else if (M(i-1,j) == min(M(i-1,j), M(i,j-1), M(i, j)) ) {
      //std::cerr << "ps I\n";
      et.push_back('I');
      i -= 1;
    }
    else if (M(i,j-1) == min(M(i-1,j), M(i, j-1), M(i, j)) ) {
      //std::cerr << "ps D\n";
      et.push_back('D');
      j -= 1;
    }
    else {
      //std::cerr << "ps U\n";
      //std::reverse(et.begin(), et.end());
      //std::cerr << "Edit transcript: \n" << et << "\n";
      // TODO: throw decent error
      exit(1);
    }
  }

  //std::cerr << "i: " << i << " j: " << j << "\n";


  /* If one string is exhausted before the other, add the necessary indels. */
  while (i > 0) {  // remaining vertical moves are insertions.
    et.push_back('I');
    i--;
  }
  while (j > 0) {  // remaining horizontal moves are deletions.
    et.push_back('D');
    j--;
  }

  //std::cerr << "Edit transcript: " << et << "\n";

  std::reverse(et.begin(), et.end());

  return {aln_score, et};
}

std::string align(const pvt::Itn &i_itn, const pvt::Itn &j_itn, pvt::aln_level_e level) {

  struct aln_args {
    pt::idx_t i_len;
    pt::idx_t j_len;
    match_res_t (*eq)(const pvt::Itn &, pt::idx_t, const pvt::Itn &, pt::idx_t);
    aln_scores_t scores;
    pvt::aln_level_e level;
  };

  auto [il, jl, eq, s, l] = ([&]() -> aln_args {
    switch (level) {
    case pvt::aln_level_e::at:
      return aln_args{i_itn.at_count(),
                      j_itn.at_count(),
                      eq_at,
                      {0, 1, 2, 1},
                      pvt::aln_level_e::at};

    case pvt::aln_level_e::step:
      return aln_args{i_itn.step_count(),
                      j_itn.step_count(),
                      eq_step,
                      {0, 1, 4, 1},
                      pvt::aln_level_e::step};
    }
  })();


  if (il == 1 && jl == 1) {
    if (eq(i_itn, 0, j_itn, 0).is_match) {
      return "M";
    }

    return "X";
  }

  auto [_, et] = global_align(i_itn, il, j_itn, jl, s, eq);

  return et;
}

} // namespace povu::align
