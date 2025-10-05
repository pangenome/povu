#ifndef PV_TINY_HPP
#define PV_TINY_HPP

#include <string_view> // for string_view

#include "povu/common/constants.hpp"	// for constants
#include "povu/graph/pvst.hpp"		// for Tree
#include "povu/graph/spanning_tree.hpp" // for Tree
#include "povu/graph/tree_utils.hpp"	// for tree_meta
#include "povu/graph/types.hpp"		// for graph

namespace povu::tiny
{
inline constexpr std::string_view MODULE = "povu::tiny";

namespace pc = povu::constants;
namespace ptu = povu::tree_utils;
namespace pst = povu::spanning_tree;
namespace pvst = povu::pvst;
namespace pgt = povu::types::graph;

void find_tiny(const pst::Tree &st, pvst::Tree &pvst, const ptu::tree_meta &tm);

} // namespace povu::tiny
#endif // PV_TINY_HPP
