#ifndef PV_PARALLEL_HPP
#define PV_PARALLEL_HPP

#include <string_view>                   // for string_view

#include "povu/graph/pvst.hpp"           // for Tree
#include "povu/graph/spanning_tree.hpp"  // for Tree
#include "povu/graph/tree_utils.hpp"     // for tree_meta
#include "povu/common/constants.hpp"     // for constants
#include "povu/graph/types.hpp"          // for graph

namespace povu::parallel
{
inline constexpr std::string_view MODULE = "povu::parallel";

namespace pc = povu::constants;
namespace ptu = povu::tree_utils;
namespace pvst = povu::pvst;
namespace pgt = povu::types::graph;
namespace pst = povu::spanning_tree;

void find_parallel(const pst::Tree &st, pvst::Tree &ft,
		   const ptu::tree_meta &tm);

} // namespace povu::parallel
#endif // PV_PARALLEL_HPP
