#include "ita/convolutions/depth_matrix.hpp"

#include <vector>

#include <liteseq/refs.h>      // for ref_walk, ref
#include <meza/pool/joint.hpp> // for joint_pool

#include "ita/convolutions/at_matrix_no_tangle.hpp" // for no_tangle
#include "ita/variation/rov.hpp"		    // for RoV
#include "povu/common/core.hpp"			    // for pt
#include "povu/graph/bidirected.hpp"		    // for VG

namespace ita::depth_matrix
{

depth_matrix comp_depth_matrix(const bd::VG &g, const ir::RoV *rov,
			       meza::pool::joint::joint_pool<qt::u32> &dm_pool)
{
	const pt::u32 I = g.get_hap_count();
	const pt::u32 J = rov->get_vertex_count();

	depth_matrix d_mat = depth_matrix(dm_pool, I, J);

	const std::vector<pt::u32> &sorted_vertices =
		rov->get_sorted_vertices();

	// auto d_mat = meza::matrix::depth_matrix{I, J};
	// auto cpy = sorted_vertices;
	// d_mat.add_col_names(std::move(cpy));

	for (pt::u32 i{}; i < I; i++) {
		for (pt::u32 j{}; j < J; j++) {
			pt::u32 v_id = sorted_vertices[j];
			pt::u32 v_idx = g.v_id_to_idx(v_id);
			const std::vector<pt::idx_t> &ref_idxs =
				g.get_vertex_ref_idxs(v_idx, i);
			pt::u32 depth = ref_idxs.size();

			d_mat.set_value(i, j, depth);

			if (depth > 1 && (j == 0 || j == J - 1))
				d_mat.set_tangled(true);

			if (depth > d_mat.get_max_depth())
				d_mat.set_max_depth(depth);
		}
	}

	return d_mat;
}

}; // namespace ita::depth_matrix
