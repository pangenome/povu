#include <set>

#include <liteseq/refs.h> // for ref_walk, ref

#include "ita/convolutions/at_matrix.hpp" // for rov_matrix_pool
#include "ita/variation/rov.hpp"	  // for RoV
#include "povu/common/core.hpp"		  // for pt
#include "povu/graph/bidirected.hpp"	  // for VG

namespace ita::at_matrix::no_tangle
{
ita::at_matrix::rov_matrix_pool
init_depth_matrices_no_tangle(const bd::VG &g, ir::RoV &rov,
			      const std::set<pt::u32> &to_call_ref_ids);

} // namespace ita::at_matrix::no_tangle
