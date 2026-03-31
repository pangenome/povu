#include <iostream>
#include <stdexcept>
#include <string>

#include <quilt/types.hpp>

#include "meza/pool/split.hpp"	  // for matrix_pool, ov_mat_t
#include "meza/shared/shared.hpp" // for layout

#include "./shared.hpp"

namespace hap_xors
{
using namespace shared;

void run(meza::pool::split::matrix_pool<qt::u8> &ov_pool)
{
	qt::u32 I = 7;
	qt::u32 J = 3;
	qt::u32 pool_offset = 0;

	std::vector<meza::pool::split::rov_mat_set> mat_sets;

	qt::u32 ctr = 0;
	qt::u32 max = 5;
	for (; ctr < max; ctr++) {
		std::cerr << "Iteration " << ctr << "\n";

		// q.pop();

		if (!ov_pool.can_allocate(
			    I, J, meza::shared::layout::DenseRowMajor)) {
			throw std::runtime_error(
				std::string("Not enough space in pool for new "
					    "marices. "
					    "Consider increasing pool size or "
					    "freeing up "
					    "space by running computations."));
		}

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

		fill_random(filter_mat);
		fill_ref_row(filter_mat, ref_mat);

		mat_sets.push_back({ref_mat, filter_mat, xor_mat, pool_offset});

		pool_offset += I * J;
	}

	ov_pool.copy_to_device();

	ov_pool.cuda_free_haps_xor();
}
}; // namespace hap_xors
