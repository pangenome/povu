#ifndef PV_MISC_HPP
#define PV_MISC_HPP

#include <string_view> // for string_view

#include "oza/graph/pvst.hpp"	       // for Tree
#include "oza/graph/spanning_tree.hpp" // for Tree
#include "oza/graph/tree_utils.hpp"    // for tree_meta

namespace oza::midi
{
inline constexpr std::string_view MODULE = "povu::misc";

namespace ptu = oza::tree_utils;

void find_midi(const pst::Tree &st, pvst::Tree &pvst, const ptu::tree_meta &tm);
} // namespace oza::midi

#endif // PV_MISC_HPP
