#include "ita/convolutions/at_matrix.hpp"

#include <thread>
#include <vector>

#include <liteseq/refs.h> // for ref_walk, ref

#include "ita/genomics/allele.hpp" // for trek
#include "ita/variation/rov.hpp"   // for RoV
#include "meza/matrix/matrix.hpp"  // for matrix2d
#include "povu/common/core.hpp"
#include "povu/graph/bidirected.hpp" // for VG
#include "povu/graph/types.hpp"	     // for ptg: or_e

namespace ita::at_matrix
{
namespace lq = liteseq;

void fill_filter_matrix(const bd::VG &g, ir::RoV &rov, pt::u32 I, pt::u32 J,
			meza::matrix::depth_matrix &m)
{
	// fill column-wise
	// row wise is more cache friendly, but needs
	// a method in bd::VG to get haplotypes per haplotype
	for (pt::u32 h_idx{}; h_idx < I; h_idx++) {
		for (pt::u32 j{}; j < J; j++) {
			pt::u32 v_id = rov.get_sorted_vertex(j);

			pt::u32 v_idx = g.v_id_to_idx(v_id);
			const std::vector<pt::idx_t> &ref_idxs =
				g.get_vertex_ref_idxs(v_idx, h_idx);
			pt::u32 depth = ref_idxs.size();

			if (depth == 0)
				continue;

			if (depth == 1) {
				const lq::ref_walk *h_w =
					g.get_ref_vec(h_idx)
						->walk;	 // the hap walk
				pt::u32 k = ref_idxs[0]; // index in the
							 // hap walk
				ptg::or_e orn =
					h_w->strands[k] ==
							lq::strand::STRAND_FWD
						? ptg::or_e::forward
						: ptg::or_e::reverse;

				switch (orn) {
				case ptg::or_e::forward:
					m.set_value(h_idx,
						    rov.get_sorted_pos(v_id),
						    1);
					break;
				case ptg::or_e::reverse:
					m.set_value(h_idx,
						    rov.get_sorted_pos(v_id),
						    2);
					break;
				}
				continue;
			}

			if (depth > 1 && (j == 0 || j == J - 1))
				m.set_tangled(true);

			if (depth > m.get_max_depth())
				m.set_max_depth(depth);

			m.set_value(h_idx, j, depth);
		}
	}
}

std::tuple<std::vector<pt::u32>, bool, pt::u32>
compute_ref_row(const bd::VG &g, ir::RoV &rov, pt::u32 ref_h_idx, pt::u32 J)
{
	std::vector<pt::u32> ref_row(J, 0);
	bool is_tangled = false;
	pt::u32 max_depth = 0;

	for (pt::u32 j{}; j < J; j++) {
		pt::u32 v_id = rov.get_sorted_vertex(j);

		pt::u32 sorted_j = rov.get_sorted_pos(v_id);

		pt::u32 v_idx = g.v_id_to_idx(v_id);
		const std::vector<pt::idx_t> &ref_idxs =
			g.get_vertex_ref_idxs(v_idx, ref_h_idx);
		pt::u32 depth = ref_idxs.size();

		if (depth == 0)
			continue;

		if (depth == 1) {
			const lq::ref_walk *h_w =
				g.get_ref_vec(ref_h_idx)->walk; // the hap walk
			pt::u32 k = ref_idxs[0];		// index in the
								// hap walk
			ptg::or_e orn =
				h_w->strands[k] == lq::strand::STRAND_FWD
					? ptg::or_e::forward
					: ptg::or_e::reverse;

			switch (orn) {
			case ptg::or_e::forward:
				ref_row[sorted_j] = 1;
				break;
			case ptg::or_e::reverse:
				ref_row[sorted_j] = 2;
				break;
			}
			continue;
		}

		if (depth > 1 && (j == 0 || j == J - 1))
			is_tangled = true;

		if (depth > max_depth)
			max_depth = depth;

		ref_row[sorted_j] = depth;
	}

	return {ref_row, is_tangled, max_depth};
}

void fill_ref_matrix(const bd::VG &g, ir::RoV &rov, pt::u32 ref_h_idx,
		     pt::u32 I, pt::u32 J, meza::matrix::depth_matrix &m)
{
	// fill column-wise
	// row wise is more cache friendly, but needs
	// a method in bd::VG to get haplotypes per haplotype

	auto [ref_row, is_tangled, max_depth] =
		compute_ref_row(g, rov, ref_h_idx, J);

	for (pt::u32 i{}; i < I; i++) {
		for (pt::u32 j{}; j < J; j++) {
			if (ref_row[j] == 0)
				continue;

			pt::u32 value = ref_row[j];
			m.set_value(i, j, value);
		}
	}

	m.set_tangled(is_tangled);
	m.set_max_depth(max_depth);
}

template <typename V>
std::vector<V> comp_row_names_fixed(pt::u32 I, V fixed_value)
{
	std::vector<V> row_names(I, fixed_value);

	return row_names;
};

std::vector<pt::u32> comp_row_names(pt::u32 I)
{
	std::vector<pt::u32> row_names;
	for (pt::u32 i{}; i < I; i++)
		row_names.push_back(i);

	return row_names;
};

matrix_pool gen_matrices(const bd::VG &g, ir::RoV &rov,
			 const std::set<pt::u32> &to_call_ref_ids, pt::u32 I,
			 pt::u32 J, std::vector<pt::u32> sorted_vertices)
{
	auto blank_matrix = meza::matrix::depth_matrix::create_full(I, J);
	blank_matrix.add_col_names(std::move(sorted_vertices));

	// duplicate the ref matrix
	auto result_matrix = blank_matrix;
	auto filter_matrix = blank_matrix;

	std::thread filter_thread(
		[&]()
		{
			fill_filter_matrix(std::cref(g), std::ref(rov), I, J,
					   std::ref(filter_matrix));
		});

	std::map<pt::u32, meza::matrix::depth_matrix> ref_matrices;
	auto ref_matrix = meza::matrix::depth_matrix::create_repeated_row(I, J);

	if (to_call_ref_ids.size() == 1)
		ref_matrices.emplace(*to_call_ref_ids.begin(), ref_matrix);
	else if (to_call_ref_ids.size() > 1)
		for (pt::u32 ref_h_idx : to_call_ref_ids)
			ref_matrices.emplace(ref_h_idx, ref_matrix);

	for (auto &[ref_h_idx, ref_matrix] : ref_matrices)
		fill_ref_matrix(g, rov, ref_h_idx, I, J, ref_matrix);

	filter_thread.join();

	return {ref_matrices, filter_matrix, result_matrix,
		filter_matrix.is_tangled()};
}

matrix_pool init_depth_matrices(const bd::VG &g, ir::RoV &rov,
				const std::set<pt::u32> &to_call_ref_ids)
{
	pt::u32 I = g.get_hap_count();	    // rows
	pt::u32 J = rov.get_vertex_count(); // cols

	const std::vector<pt::u32> &sorted_vertices = rov.get_sorted_vertices();

	matrix_pool mp =
		gen_matrices(g, rov, to_call_ref_ids, I, J, sorted_vertices);

	// auto [ref_matrices, filter_matrix, result_matrix, is_tangled] =
	//	gen_matrices(g, rov, to_call_ref_ids, I, J, sorted_vertices);

	mp.filter_matrix.add_row_names(comp_row_names(I));
	// filter_matrix.add_row_names(comp_row_names(I));

	for (pt::u32 ref_h_idx : to_call_ref_ids) {
		meza::matrix::depth_matrix &ref_matrix =
			mp.ref_matrices.at(ref_h_idx);
		std::vector<pt::u32> row_names =
			comp_row_names_fixed(I, ref_h_idx);
		ref_matrix.add_row_names(std::move(row_names));
	}

	return mp;
}

} // namespace ita::at_matrix
