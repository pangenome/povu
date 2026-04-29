#include <utility>

#include <liteseq/refs.h>     // for ref_walk, ref
#include <meza/pool/pool.hpp> // for pool
#include <quilt/types.hpp>    // for qt

#include "ita/traversals/at_matrix.hpp" // for rov_job_batch, mat3, mat3_item, hap2loop
#include "ita/traversals/depth_matrix.hpp" // for depth_matrix, comp_depth_matrix
#include "ita/variation/rov.hpp"	   // for RoV

namespace ita::at_matrix::no_tangle
{
namespace lq = liteseq;
using meza::matrix::depth_matrix;

constexpr meza::pool::pool_region region_ref =
	meza::pool::pool_region::Reference;

constexpr meza::pool::pool_region region_filter =
	meza::pool::pool_region::Filter;

constexpr meza::pool::pool_region region_xor = meza::pool::pool_region::Xor;

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
		    const std::set<qt::u32> &to_call_ref_ids,
		    const ita::depth_matrix::depth_matrix dm, pool_t &p,
		    rov_job_batch &batch)
{
	qt::u32 I = dm.rows();
	qt::u32 J = dm.cols();

	rov_job j{rov, I, {}};

	for (const qt::u32 ref_h_idx : to_call_ref_ids) {
		auto ref_mat = p.alloc_ov_matrix<std::string, std::string>(
			I, J, region_ref);

		auto filter_mat = p.alloc_ov_matrix<std::string, std::string>(
			I, J, region_filter);

		auto xor_mat = p.alloc_ov_matrix<std::string, std::string>(
			I, J, region_xor);

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
