#ifndef OZ_CONCEALED_HPP
#define OZ_CONCEALED_HPP

#include <string_view> // for string_view

#include <quilt/constants.hpp>	 // for
#include <quilt/graph_types.hpp> // for v_end_e, side_n_id_t, side_n_idx_t

#include "oza/graph/pvst.hpp"	       // for Tree
#include "oza/graph/spanning_tree.hpp" // for Tree
#include "oza/graph/tree_utils.hpp"    // for tree_meta

namespace oza::concealed
{
inline constexpr std::string_view MODULE = "povu::concealed";

namespace pst = oza::spanning_tree;
namespace pvst = oza::pvst;
namespace ptu = oza::tree_utils;

void find_concealed(const pst::Tree &st, pvst::Tree &ft,
		    const ptu::tree_meta &tm);
} // namespace oza::concealed
#endif // OZ_CONCEALED_HPP
