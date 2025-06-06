#ifndef PV_FLUBBLES_HPP
#define PV_FLUBBLES_HPP

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <deque>
#include <iostream>
#include <limits>
#include <map>
#include <stack>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "../../src/cli/app.hpp"
#include "../common/tree_utils.hpp"
#include "../common/types.hpp"
#include "../graph/spanning_tree.hpp"
#include "../graph/tree.hpp"

namespace povu::flubbles {

#define MODULE "povu::graph::flubble_tree"

namespace pgt = povu::graph_types;
using namespace povu::graph_types;

namespace pt = povu::types;
namespace pvtr = povu::tree;
namespace pc = povu::constants;
namespace pvst = povu::types::pvst;
namespace pst = povu::spanning_tree;
namespace ptu = povu::tree_utils;

// a hairpin boundary
struct boundary {
  std::size_t b1;
  std::size_t b2;
};

const pt::idx_t EXPECTED_HAIRPIN_COUNT{1000}; // Expected number of hairpins
const boundary NULL_BOUNDARY{.b1 = pc::INVALID_IDX, .b2 = pc::INVALID_IDX};

// orientation, id, class
struct oic_t {
  pgt::or_e orientation;
  pt::id_t id;
  pt::id_t st_idx;
  pt::idx_t cls;
};

struct eq_class_stack_t {
  std::vector<oic_t> s; // stack of equivalence classes
  std::vector<pt::idx_t>
      next_seen; // next seen index for each equivalence class

  // constructor
  eq_class_stack_t(pt::idx_t exp_size) {
    s.reserve(exp_size);         // reserve some space for the stack
    next_seen.reserve(exp_size); // reserve some space for the next seen vector
  }
};

/**
 * @brief Generate flubble tree from spanning tree
 *
 */
pvtr::Tree<pvst::Vertex> find_flubbles(pst::Tree &t, const core::config &app_config);
} // namespace povu::flubbles
#endif
