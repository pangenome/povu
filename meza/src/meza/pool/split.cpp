#include "meza/pool/split.hpp"

#include <driver_types.h>
#include <set>

// #if MEZA_USE_CUDA
// #include <cuda_runtime.h>
// #include <cuda_runtime_api.h>
// #endif

#include "quilt/types.hpp"

namespace meza::pool::split
{

haps_comp_set handle_set(matrix_pool<qt::u8> &ov_pool,
			 const ov_mat_t &filter_mat, qt::u32 pool_offset)
{
	hap_comp_matrix comp_mat(ov_pool.get_haps_xor(), ov_pool.get_haps_sum(),
				 filter_mat, pool_offset);

	// TODO: parallelise the run_in_haps calls
	comp_mat.run_in_haps(ov_pool, comparison_op::XOR);
	comp_mat.run_in_haps(ov_pool, comparison_op::SUM);

	std::set<qt::up_t<qt::u32>> reversals = comp_mat.find_reversals();

	auto [matches, mismatches] = comp_mat.explore_pairs();

	return {reversals, matches, mismatches};
}

} // namespace meza::pool::split
