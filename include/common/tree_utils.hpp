#ifndef PV_TREE_UTILS_HPP
#define PV_TREE_UTILS_HPP

#include <algorithm>
#include <cassert>
#include <list>
#include <map>
#include <set>
#include <stack>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../graph/spanning_tree.hpp"
//#include "../graph/tree.hpp"
#include "./types.hpp"
namespace povu::tree_utils {

#define MODULE "povu::tree_utils"

namespace pgt = povu::graph_types;
using namespace povu::graph_types;

namespace pt = povu::types;
  //namespace pvtr = povu::tree;
namespace pc = povu::constants;
  //namespace pvst = povu::types::pvst;
namespace pst = povu::spanning_tree;

// key is edge idx of an edge and value is the v idx of the next braching vtx or
// of a leaf
using EdgeToBranch = std::map<pt::idx_t, pt::idx_t>;
struct BranchingMeta {
  // the black edge idx
  pt::idx_t black_e_idx = pc::INVALID_IDX;

  // branches sorted from the highest to lowest but if any child is black then
  // add first
  std::vector<pt::idx_t> sorted_br;

  // the vector should be sorted ...
  EdgeToBranch edge_data;
};

// key is a vertex idx and value is the branching meta data
using BranchDesc = std::map<pt::idx_t, BranchingMeta>;

BranchDesc br_desc(const pst::Tree &st);
} // namespace povu::tree_utils

#endif
