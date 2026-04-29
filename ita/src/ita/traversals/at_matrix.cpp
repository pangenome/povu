#include "ita/traversals/at_matrix.hpp"

#include <meza/pool/pool.hpp> // for matrix_pool

#include <liteseq/refs.h>  // for ref_walk, ref
#include <quilt/types.hpp> // for qt

#include "ita/traversals/at_matrix_no_tangle.hpp" // for no_tangle, from_no_tangle
#include "ita/traversals/at_matrix_tangled.hpp" // for init_tangled_depth_matrices
#include "ita/traversals/depth_matrix.hpp" // for depth_matrix, comp_depth_matrix
#include "ita/traversals/untangle.hpp"	   // for untangle
#include "ita/variation/rov.hpp"	   // for RoV
#include "povu/graph/bidirected.hpp"	   // for VG

namespace ita::at_matrix
{
using pool_t = meza::pool::pool<qt::u8, qt::u32>;

void init_pool(const bd::VG &g, const ir::RoV *rov,
	       const std::set<qt::u32> &to_call_ref_ids,
	       const ita::depth_matrix::depth_matrix &dm, pool_t &p,
	       rov_job_batch &batch)
{
	if (dm.is_tangled()) {
		ita::traversals::untangle::aln_chain c =
			ita::traversals::untangle::untangle(g, to_call_ref_ids,
							    *rov);
		ita::at_matrix::tangled::from_tangled(g, rov, to_call_ref_ids,
						      p, c, batch);
	}
	else {
		ita::at_matrix::no_tangle::from_no_tangle(rov, to_call_ref_ids,
							  dm, p, batch);
	}
}

} // namespace ita::at_matrix
