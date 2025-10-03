#ifndef PV_FLUBBLES_HPP
#define PV_FLUBBLES_HPP

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <deque>
#include <iostream>
#include <limits>
#include <map>
#include <stack>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

"#include "povu/common/app.hpp"
"#include "povu/common/compat.hpp"
"#include "povu/graph/pvst.hpp"
"#include "povu/graph/spanning_tree.hpp"
"#include "povu/graph/tree_utils.hpp"

namespace povu::flubbles
{
inline constexpr std::string_view MODULE = "povu::graph::flubble_tree";

namespace pc = povu::constants;
namespace ptu = povu::tree_utils;
namespace pvst = povu::pvst;
namespace pgt = povu::types::graph;
namespace pst = povu::spanning_tree;

// a hairpin boundary
struct boundary {
	std::size_t b1;
	std::size_t b2;
};

const pt::idx_t EXPECTED_HAIRPIN_COUNT{1000}; // Expected number of hairpins
const boundary NULL_BOUNDARY{pc::INVALID_IDX, pc::INVALID_IDX};

// orientation, id, class
struct oic_t {
	pgt::or_e orientation;
	pt::id_t id;
	pt::id_t st_idx;
	pt::idx_t cls;
};

struct eq_class_stack_t {
	std::vector<oic_t> s; // stack of equivalence classes
	std::vector<pt::idx_t>
		next_seen; // next seen index for each equivalence class

	// constructor
	eq_class_stack_t(pt::idx_t exp_size)
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
} // namespace povu::flubbles
#endif
