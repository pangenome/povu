#ifndef PV_SMOTHERED_HPP
#define PV_SMOTHERED_HPP

#include <string_view> // for string_view

#include "povu/common/constants.hpp"	// for constants
#include "povu/graph/pvst.hpp"		// for Tree
#include "povu/graph/spanning_tree.hpp" // for Tree
#include "povu/graph/tree_utils.hpp"	// for tree_meta
#include "povu/graph/types.hpp"		// for graph

namespace povu::smothered
{
inline constexpr std::string_view MODULE = "povu::smothered";

namespace pgt = povu::types::graph;
namespace pc = povu::constants;
namespace pst = povu::spanning_tree;
namespace pvst = povu::pvst;
namespace pc = povu::constants;
namespace ptu = povu::tree_utils;

void find_smothered(const pst::Tree &st, pvst::Tree &ft,
		    const ptu::tree_meta &tm);
} // namespace povu::smothered

#endif // PV_SMOTHERED_HPP
