#ifndef OZ_SMOTHERED_HPP
#define OZ_SMOTHERED_HPP

#include <string_view> // for string_view

#include "oza/graph/pvst.hpp"	       // for Tree
#include "oza/graph/spanning_tree.hpp" // for Tree
#include "oza/graph/tree_utils.hpp"    // for tree_meta

namespace oza::smothered
{

namespace ptu = oza::tree_utils;

void find_smothered(const pst::Tree &st, pvst::Tree &ft,
		    const ptu::tree_meta &tm);
} // namespace oza::smothered

#endif // OZ_SMOTHERED_HPP
