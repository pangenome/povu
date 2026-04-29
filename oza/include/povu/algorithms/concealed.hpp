#ifndef PV_CONCEALED_HPP
#define PV_CONCEALED_HPP

#include <string_view> // for string_view

#include "povu/common/constants.hpp"	// for constants
#include "povu/graph/pvst.hpp"		// for Tree
#include "povu/graph/spanning_tree.hpp" // for Tree
#include "povu/graph/tree_utils.hpp"	// for tree_meta
#include "povu/graph/types.hpp"		// for graph

namespace oza::concealed
{

inline constexpr std::string_view MODULE = "povu::concealed";

namespace pgt = povu::types::graph;
namespace pc = povu::constants;
namespace pst = oza::spanning_tree;
namespace pvst = oza::pvst;
namespace pc = povu::constants;
namespace ptu = oza::tree_utils;

void find_concealed(const pst::Tree &st, pvst::Tree &ft,
		    const ptu::tree_meta &tm);
} // namespace oza::concealed
#endif // PV_CONCEALED_HPP
