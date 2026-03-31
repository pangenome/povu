#include <gtest/gtest.h>

#include <meza/pool/split.hpp>
#include <quilt/types.hpp>

TEST(SharedPool, InitPool)
{
	auto &ov_pool = meza::pool::split::matrix_pool<qt::u8>::init();

	// fixed matrix size
	qt::u32 I = 5;
	qt::u32 J = 3;

	qt::u32 max_matrices = 5;
	qt::u32 ctr = 0;

	while (ctr < max_matrices) {
		if (!ov_pool.can_allocate(
			    I, J, meza::pool::split::layout::DenseRowMajor)) {
			ov_pool.reset();
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

		ctr++;
	}

	ASSERT_FALSE(ov_pool.empty());
}
