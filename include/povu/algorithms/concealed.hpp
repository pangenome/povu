#ifndef PV_CONCEALED_HPP
#define PV_CONCEALED_HPP

#include <string_view> // for string_view

#include "povu/common/constants.hpp"	// for constants
#include "povu/graph/pvst.hpp"		// for Tree
#include "povu/graph/spanning_tree.hpp" // for Tree
#include "povu/graph/tree_utils.hpp"	// for tree_meta
#include "povu/graph/types.hpp"		// for graph

namespace povu::concealed
{

inline constexpr std::string_view MODULE = "povu::concealed";

namespace pgt = povu::types::graph;
namespace pc = povu::constants;
namespace pst = povu::spanning_tree;
namespace pvst = povu::pvst;
namespace pc = povu::constants;
namespace ptu = povu::tree_utils;

void find_concealed(const pst::Tree &st, pvst::Tree &ft,
		    const ptu::tree_meta &tm);
} // namespace povu::concealed
#endif // PV_CONCEALED_HPP
