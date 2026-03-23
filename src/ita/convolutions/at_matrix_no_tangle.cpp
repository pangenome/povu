#include <vector>

#include <liteseq/refs.h>      // for ref_walk, ref
#include <meza/pool/split.hpp> // for matrix_pool
#include <quilt/types.hpp>

#include "ita/convolutions/at_matrix.hpp"    //
#include "ita/convolutions/depth_matrix.hpp" // for depth_matrix, comp_depth_matrix
#include "ita/variation/rov.hpp"	     // for RoV
#include "povu/common/core.hpp"		     // for pt
#include "povu/graph/bidirected.hpp"	     // for VG
#include "povu/graph/types.hpp"		     // for ptg: or_e

namespace ita::at_matrix::no_tangle
{
namespace lq = liteseq;
using meza::matrix::depth_matrix;

void fill_filter_matrix(const bd::VG &g, const ir::RoV &rov, pt::u32 I,
			pt::u32 J, meza::matrix::depth_matrix &m)
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
compute_ref_row(const bd::VG &g, const ir::RoV &rov, pt::u32 ref_h_idx,
		pt::u32 J)
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

void fill_ref_matrix(const bd::VG &g, const ir::RoV &rov, pt::u32 ref_h_idx,
		     pt::u32 I, pt::u32 J, meza::matrix::ref_matrix &m)
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

// matrix_pool gen_matrices(const bd::VG &g, const ir::RoV &rov,
//			 const std::set<pt::u32> &to_call_ref_ids, pt::u32 I,
//			 pt::u32 J, std::vector<pt::u32> sorted_vertices)
// {
//	auto blank_matrix = meza::matrix::depth_matrix(I, J);
//	blank_matrix.add_col_names(std::move(sorted_vertices));

//	auto result_matrix = blank_matrix;
//	auto filter_matrix = blank_matrix;

//	std::thread filter_thread(
//		[&]()
//		{
//			fill_filter_matrix(std::cref(g), std::ref(rov), I, J,
//					   std::ref(filter_matrix));
//		});

//	std::map<pt::u32, meza::matrix::ref_matrix> ref_matrices;
//	auto ref_matrix = meza::matrix::ref_matrix{I, J};

//	if (to_call_ref_ids.size() == 1)
//		ref_matrices.emplace(*to_call_ref_ids.begin(), ref_matrix);
//	else if (to_call_ref_ids.size() > 1)
//		for (pt::u32 ref_h_idx : to_call_ref_ids)
//			ref_matrices.emplace(ref_h_idx, ref_matrix);

//	for (auto &[ref_h_idx, ref_matrix] : ref_matrices)
//		fill_ref_matrix(g, rov, ref_h_idx, I, J, ref_matrix);

//	filter_thread.join();

//	return {ref_matrices, filter_matrix, result_matrix, {}, 0, false, I, J};
// }

// rov_matrix_pool
// init_depth_matrices_no_tangle(const bd::VG &g, const ir::RoV *rov,
//			      const std::set<pt::u32> &to_call_ref_ids)
// {
//	rov_matrix_pool rov_mp{*rov};
//	pt::u32 I = g.get_hap_count();	     // rows
//	pt::u32 J = rov->get_vertex_count(); // cols

//	const std::vector<pt::u32> &sorted_vertices =
//		rov->get_sorted_vertices();

//	matrix_pool mp =
//		gen_matrices(g, *rov, to_call_ref_ids, I, J, sorted_vertices);

//	mp.I = I;
//	mp.J = J;

//	mp.filter_matrix.add_row_names(comp_row_names(I));

//	for (pt::u32 ref_h_idx : to_call_ref_ids) {
//		meza::matrix::ref_matrix &ref_matrix =
//			mp.ref_matrices.at(ref_h_idx);

//		std::vector<pt::u32> row_names =
//			comp_row_names_fixed(I, ref_h_idx);

//		ref_matrix.add_row_names(std::move(row_names));
//	}

//	rov_mp.pools.emplace_back(std::move(mp));

//	return rov_mp;
// }

void populate_filter(const qt::u32 I, const qt::u32 J,
		     const ita::depth_matrix::depth_matrix &dm,
		     ov_mat_t &filter_mat)
{
	for (qt::u32 i{}; i < I; i++) {
		for (qt::u32 j{}; j < J; j++) {
			qt::u8 val = dm.get_value(i, j);
			filter_mat.set_value(i, j, val);
		}
	}
}

void populate_ref(const qt::u32 I, const qt::u32 J, qt::u32 ref_i,
		  const ita::depth_matrix::depth_matrix &dm, ov_mat_t &ref_mat)
{
	for (qt::u32 i{}; i < I; i++) {
		for (qt::u32 j{}; j < J; j++) {
			qt::u8 val = dm.get_value(ref_i, j);
			ref_mat.set_value(i, j, val);
		}
	}
}

void from_no_tangle(const ir::RoV *rov,
		    const std::set<pt::u32> &to_call_ref_ids,
		    const ita::depth_matrix::depth_matrix dm,
		    meza::pool::split::matrix_pool<qt::u8> &ov_pool,
		    rov_job_batch &batch)
{
	qt::u32 I = dm.rows();
	qt::u32 J = dm.cols();

	rov_job j{rov, {}};

	for (const pt::u32 ref_h_idx : to_call_ref_ids) {
		auto ref_mat =
			ov_pool.alloc_ov_matrix<std::string, std::string>(
				I, J,
				meza::pool::split::pool_region::Reference);

		auto filter_mat =
			ov_pool.alloc_ov_matrix<std::string, std::string>(
				I, J, meza::pool::split::pool_region::Filter);

		auto xor_mat =
			ov_pool.alloc_ov_matrix<std::string, std::string>(
				I, J, meza::pool::split::pool_region::Xor);

		qt::u32 filter_size = filter_mat.base().size();

		populate_filter(I, J, dm, filter_mat);
		populate_ref(I, J, ref_h_idx, dm, ref_mat);

		mat3 m{std::move(ref_mat),
		       std::move(filter_mat),
		       std::move(xor_mat),
		       batch.pool_j_offset,
		       I,
		       J};

		mat3_item item{std::move(m)};

		j.add_item(ref_h_idx, std::move(item));

		batch.pool_j_offset += filter_size;
	}

	batch.add(std::move(j));
}

} // namespace ita::at_matrix::no_tangle
