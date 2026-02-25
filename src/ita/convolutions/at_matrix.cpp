#include "ita/convolutions/at_matrix.hpp"

#include <thread>
#include <vector>

#include <liteseq/refs.h> // for ref_walk, ref

#include "ita/convolutions/at_matrix_no_tangle.hpp" // for no_tangle
#include "ita/convolutions/at_matrix_tangled.hpp" // for init_tangled_depth_matrices
#include "ita/traversals/untangle.hpp"		  // for untangle
#include "ita/variation/rov.hpp"		  // for RoV
#include "meza/matrix/matrix.hpp"		  // for matrix2d
#include "povu/common/core.hpp"			  // for pt
#include "povu/graph/bidirected.hpp"		  // for VG

namespace ita::at_matrix
{
using ita::at_matrix::no_tangle::init_depth_matrices_no_tangle;

meza::matrix::depth_matrix
comp_depth_matrix(const bd::VG &g, const ir::RoV &rov,
		  const std::vector<pt::u32> &sorted_vertices)
{
	const pt::u32 I = g.get_hap_count();
	const pt::u32 J = rov.get_vertex_count();

	auto d_mat = meza::matrix::depth_matrix{I, J};
	auto cpy = sorted_vertices;
	d_mat.add_col_names(std::move(cpy));

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

	// std::cerr << rov.as_str() << " Depth matrix\n";
	// d_mat.base().dbg_print(std::cerr);

	return d_mat;
}

rov_matrix_pool init_tangled(const bd::VG &g, ir::RoV &rov,
			     const std::set<pt::u32> &to_call_ref_ids)
{
	ita::traversals::untangle::aln_chain c =
		ita::traversals::untangle::untangle(g, to_call_ref_ids, rov);
	ita::at_matrix::rov_matrix_pool tangled_rov_mp =
		ita::at_matrix::tangled::init_tangled_depth_matrices(
			g, rov, to_call_ref_ids, c);
	return tangled_rov_mp;
}

inline rov_matrix_pool
init_not_tangled(const bd::VG &g, ir::RoV &rov,
		 const std::set<pt::u32> &to_call_ref_ids)
{
	return init_depth_matrices_no_tangle(g, rov, to_call_ref_ids);
}

rov_matrix_pool init_depth_matrices(const bd::VG &g, ir::RoV &rov,
				    const std::set<pt::u32> &to_call_ref_ids)
{
	const std::vector<pt::u32> &sorted_vertices = rov.get_sorted_vertices();
	meza::matrix::depth_matrix d_mat =
		comp_depth_matrix(g, rov, sorted_vertices);
	return d_mat.is_tangled() ? init_tangled(g, rov, to_call_ref_ids)
				  : init_not_tangled(g, rov, to_call_ref_ids);
}

} // namespace ita::at_matrix
