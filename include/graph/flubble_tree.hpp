#ifndef FLUBBLE_TREE_HPP
#define FLUBBLE_TREE_HPP


#include <vector>

#include "../graph/spanning_tree.hpp"
#include "../graph/tree.hpp"
#include "../common/types.hpp"

// rename this namespace to PVST
namespace povu::graph::flubble_tree {

#define MODULE "povu::graph::flubble_tree"

namespace pvtr = povu::tree;
namespace pst = povu::spanning_tree;
namespace pgt = povu::graph_types;
namespace pvst = povu::types::pvst;

/**
 * @brief
 */
std::vector<graph_types::canonical_sese> find_seses(pst::Tree& t);

/**
 * @brief Generate flubble tree from spanning tree
 *
 */
pvtr::Tree<pvst::Vertex> st_to_ft(pst::Tree &t);
}


#endif
