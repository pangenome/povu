#ifndef PV_PARALLEL_HPP
#define PV_PARALLEL_HPP

#include <string_view> // for string_view

#include "oza/graph/pvst.hpp"	       // for Tree
#include "oza/graph/spanning_tree.hpp" // for Tree
#include "oza/graph/tree_utils.hpp"    // for tree_meta

namespace oza::parallel
{
inline constexpr std::string_view MODULE = "povu::parallel";

namespace ptu = oza::tree_utils;

void find_parallel(const pst::Tree &st, pvst::Tree &ft,
		   const ptu::tree_meta &tm);

} // namespace oza::parallel
#endif // PV_PARALLEL_HPP
