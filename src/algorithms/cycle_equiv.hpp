#ifndef VST_HPP
#define VST_HPP

#include <cstddef>

#include "../core/constants.hpp"
#include "../graph/spanning_tree.hpp"

namespace algorithms {
void cycle_equiv(spanning_tree::Tree &t);

struct branch {
  std::size_t id {core::constants::UNDEFINED_SIZE_T};
  std::vector<std::size_t> data;
  bool is_nesting {false};
};

/**
 * A struct to hold the branch vectors of a vertex
 */
struct branch_vecs {
  std::vector<branch> branches;
  std::vector<std::size_t> merged;
  std::size_t parent_id {core::constants::UNDEFINED_SIZE_T};

  // std::size_t nesting_branch;

  // std::vector<std::size_t> left;
  // std::vector<std::size_t> right;

  //std::vector<std::size_t> children;
  //std::size_t left_id {core::constants::UNDEFINED_SIZE_T};
  //std::size_t right_id {core::constants::UNDEFINED_SIZE_T};
};

std::vector<std::size_t> compute_eq_class_stack(spanning_tree::Tree st);


} // namespace algorithms



#endif
