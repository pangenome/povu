#include "ita/traversals/depth_matrix.hpp"

#include <vector>

#include <liteseq/refs.h>	    // for ref_walk, ref
#include <meza/pool/joint.hpp>	    // for joint_pool
#include <meza/pool/pool.hpp>	    // for pool
#include <oza/graph/bidirected.hpp> // for VG
#include <quilt/types.hpp>	    // for qt

// #include "ita/traversals/at_matrix_no_tangle.hpp" // for no_tangle
#include "ita/variation/rov.hpp" // for RoV

namespace ita::depth_matrix
{

depth_matrix comp_depth_matrix(const bd::VG &g, const ir::RoV *rov, pool_t &p)
{
	const qt::u32 I = g.get_hap_count();
	const qt::u32 J = rov->get_vertex_count();

	// meza::pool::joint::full_view<qt::u32> dm_view =
	//	p.alloc_depth_matrix(I, J);

	depth_matrix d_mat = depth_matrix{p.alloc_depth_matrix(I, J)};

	const std::vector<qt::u32> &sorted_vertices =
		rov->get_sorted_vertices();

	// auto d_mat = meza::matrix::depth_matrix{I, J};
	// auto cpy = sorted_vertices;
	// d_mat.add_col_names(std::move(cpy));

	for (qt::u32 i{}; i < I; i++) {
		for (qt::u32 j{}; j < J; j++) {
			qt::u32 v_id = sorted_vertices[j];
			qt::u32 v_idx = g.v_id_to_idx(v_id);
			const std::vector<qt::idx_t> &ref_idxs =
				g.get_vertex_ref_idxs(v_idx, i);
			qt::u32 depth = ref_idxs.size();

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
