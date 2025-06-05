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

#include "../common/types.hpp"
#include "../common/tree_utils.hpp"
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

/**
 * @brief Generate flubble tree from spanning tree
 *
 */
pvtr::Tree<pvst::Vertex> find_flubbles(const pst::Tree &t);
} // namespace povu::flubbles
#endif
