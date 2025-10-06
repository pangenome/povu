#ifndef PV_MISC_HPP
#define PV_MISC_HPP

#include <string_view> // for string_view

#include "povu/common/constants.hpp"	// for constants
#include "povu/graph/pvst.hpp"		// for Tree
#include "povu/graph/spanning_tree.hpp" // for Tree
#include "povu/graph/tree_utils.hpp"	// for tree_meta
#include "povu/graph/types.hpp"		// for graph

namespace povu::midi
{
inline constexpr std::string_view MODULE = "povu::misc";

namespace pgt = povu::types::graph;
namespace pc = povu::constants;
namespace pst = povu::spanning_tree;
namespace pvst = povu::pvst;
namespace pc = povu::constants;
namespace ptu = povu::tree_utils;

void find_midi(const pst::Tree &st, pvst::Tree &pvst, const ptu::tree_meta &tm);
} // namespace povu::midi

#endif // PV_MISC_HPP
