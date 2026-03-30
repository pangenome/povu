#include <utility>

#include <liteseq/refs.h>      // for ref_walk, ref
#include <meza/pool/split.hpp> // for matrix_pool
#include <quilt/types.hpp>

#include "ita/traversals/at_matrix.hpp"	   //
#include "ita/traversals/depth_matrix.hpp" // for depth_matrix, comp_depth_matrix
#include "ita/variation/rov.hpp"	   // for RoV
#include "povu/common/core.hpp"		   // for pt

namespace ita::at_matrix::no_tangle
{
namespace lq = liteseq;
using meza::matrix::depth_matrix;

void populate_filter2(const ita::depth_matrix::depth_matrix &dm,
		      ov_mat_t &filter_mat)
{
	auto [src_ptr, len] = dm.get_slice();

	filter_mat.copy_slice(src_ptr, len);
}

void populate_ref2(const qt::u32 I, const qt::u32 J, qt::u32 ref_i,
		   const ita::depth_matrix::depth_matrix &dm, ov_mat_t &ref_mat)
{
	auto [src_ptr, len] = dm.get_slice();
	qt::u32 ref_row = ref_i * J; // offset
	const qt::u32 *src_row_ptr = src_ptr + ref_row;
	qt::u32 repeat = I; // number of times to repeat the ref row

	ref_mat.copy_slice(src_row_ptr, J, repeat);
}

void from_no_tangle(const ir::RoV *rov,
		    const std::set<pt::u32> &to_call_ref_ids,
		    const ita::depth_matrix::depth_matrix dm,
		    meza::pool::split::matrix_pool<qt::u8> &ov_pool,
		    rov_job_batch &batch)
{
	qt::u32 I = dm.rows();
	qt::u32 J = dm.cols();

	rov_job j{rov, I, {}};

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

		populate_filter2(dm, filter_mat);
		populate_ref2(I, J, ref_h_idx, dm, ref_mat);

		mat3 m{std::move(ref_mat),
		       std::move(filter_mat),
		       std::move(xor_mat),
		       batch.get_pool_j_offset(),
		       I,
		       J};

		mat3_item item{std::move(m)};

		j.add_item(ref_h_idx, std::move(item));
		qt::u32 offset = batch.get_pool_j_offset() + filter_size;
		batch.set_pool_j_offset(offset);
	}

	batch.add(std::move(j));
}

} // namespace ita::at_matrix::no_tangle
