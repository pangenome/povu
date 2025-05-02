#ifndef FLUBBLE_TREE_HPP
#define FLUBBLE_TREE_HPP


#include <vector>

#include "../graph/spanning_tree.hpp"
#include "../graph/tree.hpp"
#include "../common/types.hpp"

namespace povu::graph::flubble_tree {

#define MODULE "povu::graph::flubble_tree"

namespace pvtr = povu::tree;
namespace pst = povu::spanning_tree;
namespace pgt = povu::graph_types;
  // TODO: rename to something related to reflect return type

/**
 * @brief
 */
std::vector<graph_types::canonical_sese> find_seses(pst::Tree& t);

/**
 * @brief spanning tree vertex to bubble tree
 *
 */
pvtr::Tree<pgt::flubble> st_to_ft(pst::Tree& t);
  //std::vector<pgt::flubble> enumerate(pst::Tree& t);
}


#endif
