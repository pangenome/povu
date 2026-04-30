#ifndef PV_FLUBBLES_HPP
#define PV_FLUBBLES_HPP

#include <cstddef>     // for size_t
#include <string_view> // for string_view
#include <vector>      // for vector

#include <quilt/graph_types.hpp> // for v_end_e, side_n_id_t, side_n_idx_t
#include <quilt/types.hpp>	 // for qt

#include "povu/common/app.hpp" // for config
// #include "povu/common/constants.hpp"	// for INVALID_IDX
#include "povu/graph/pvst.hpp"		// for Tree
#include "povu/graph/spanning_tree.hpp" // for Tree
#include "povu/graph/tree_utils.hpp"	// for tree_utils

namespace oza::flubbles
{
inline constexpr std::string_view MODULE = "povu::graph::flubble_tree";

namespace ptu = oza::tree_utils;
namespace pgt = quilt::types::graph;

// a hairpin boundary
struct boundary {
	std::size_t b1;
	std::size_t b2;
};

const qt::idx_t EXPECTED_HAIRPIN_COUNT{1000}; // Expected number of hairpins
const boundary NULL_BOUNDARY{pc::INVALID_IDX, pc::INVALID_IDX};

// orientation, id, class
struct oic_t {
	pgt::or_e orientation;
	qt::id_t id;
	qt::id_t st_idx;
	qt::idx_t cls;
};

struct eq_class_stack_t {
	std::vector<oic_t> s; // stack of equivalence classes
	std::vector<qt::idx_t>
		next_seen; // next seen index for each equivalence class

	// constructor
	eq_class_stack_t(qt::idx_t exp_size)
	{
		s.reserve(exp_size);	     // reserve some space for the stack
		next_seen.reserve(exp_size); // reserve some space for the next
					     // seen vector
	}
};

/**
 * @brief Generate flubble tree from spanning tree
 *
 */
pvst::Tree find_flubbles(pst::Tree &t, const core::config &app_config);
} // namespace oza::flubbles
#endif
