#include <gtest/gtest.h>

#include <meza/pool/split.hpp>
#include <quilt/types.hpp>

using qt::u32;
using qt::u8;

TEST(SplitPool, Init)
{
	auto &ov_pool = meza::pool::split::matrix_pool<u8>::init();
	ASSERT_TRUE(ov_pool.empty());
}

TEST(SplitPool, InitVariedSize)
{
	constexpr auto DRM = meza::pool::split::layout::DenseRowMajor;

	// fixed matrix size
	// random I and J values between 1 and 10 for each matrix
	// random value between 1 and 10
	u32 MAX_I{10};
	u32 MAX_J{10};
	u32 MAX_MATRICES{10};

	u32 I = std::rand() % MAX_I + 1;
	u32 J = std::rand() % MAX_J + 1;

	u32 max_matrices = std::rand() % MAX_MATRICES;

	u32 ctr{};
	u32 expected_used{};

	auto &ov_pool = meza::pool::split::matrix_pool<u8>::init();

	while (ctr < max_matrices) {
		if (!ov_pool.can_allocate(I, J, DRM))
			ov_pool.reset();

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

		auto [ref_used, filter_used, xor_used] = ov_pool.used();

		ctr++;

		expected_used += I * J;

		ASSERT_EQ(ref_used, expected_used);
		ASSERT_EQ(filter_used, expected_used);
		ASSERT_EQ(xor_used, expected_used);
	}

	ASSERT_FALSE(ov_pool.empty());

	// we don't expect it to be full when MAX_MATRICES is small
	ASSERT_FALSE(ov_pool.is_full());
}
