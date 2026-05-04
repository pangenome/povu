#ifndef OZ_TINY_HPP
#define OZ_TINY_HPP

#include <string_view> // for string_view

#include "oza/graph/pvst.hpp"	       // for Tree
#include "oza/graph/spanning_tree.hpp" // for Tree
#include "oza/graph/tree_utils.hpp"    // for tree_meta

namespace oza::tiny
{
inline constexpr std::string_view MODULE = "povu::tiny";

namespace ptu = oza::tree_utils;

void find_tiny(const pst::Tree &st, pvst::Tree &pvst, const ptu::tree_meta &tm);

} // namespace oza::tiny
#endif // OZ_TINY_HPP
