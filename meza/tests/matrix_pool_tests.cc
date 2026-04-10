#include <gtest/gtest.h>

#include <quilt/types.hpp>

#include <meza/pool/matrix_pool.hpp>
#include <meza/pool/split_pool_types.hpp>

#if MEZA_USE_CUDA
#include "meza/pool/matrix_pool_cuda.cuh"
#endif

using qt::u32;
using qt::u8;

/* === init pool tests === */

TEST(SplitPool, InitCpuOnly)
{
	auto p = meza::pool::matrix_pool<u8>::create_from_megabytes();
	ASSERT_TRUE(p.empty());
}

TEST(SplitPool, InitCuda)
{
#if MEZA_USE_CUDA
	auto p = meza::pool::matrix_pool<u8>::create_from_megabytes();
	auto p_cuda = meza::pool::matrix_pool_cuda<u8>(p);
	ASSERT_TRUE(p_cuda.base().empty());
#endif
}

/* === alloc single matrix tests === */

TEST(SplitPool, AllocCpuSingleMatrix)
{
	auto p = meza::pool::matrix_pool<u8>::create_from_megabytes();

	u32 MAX_I{10};
	u32 MAX_J{10};

	u32 I = std::rand() % MAX_I + 1;
	u32 J = std::rand() % MAX_J + 1;

	auto ref_mat = p.alloc_ov_matrix<std::string, std::string>(
		I, J, meza::pool::pool_region::Reference);

	auto filter_mat = p.alloc_ov_matrix<std::string, std::string>(
		I, J, meza::pool::pool_region::Filter);

	auto xor_mat = p.alloc_ov_matrix<std::string, std::string>(
		I, J, meza::pool::pool_region::Xor);

	auto [ref_used, filter_used, xor_used] = p.used();

	ASSERT_EQ(ref_used, I * J);
	ASSERT_EQ(filter_used, I * J);
	ASSERT_EQ(xor_used, I * J);

	ASSERT_FALSE(p.empty());
}

TEST(SplitPool, AllocCudaSingleMatrix)
{
#if MEZA_USE_CUDA
	auto p = meza::pool::matrix_pool<u8>::create_from_megabytes();
	auto p_cuda = meza::pool::matrix_pool_cuda<u8>(p);

	u32 MAX_I{10};
	u32 MAX_J{10};

	u32 I = std::rand() % MAX_I + 1;
	u32 J = std::rand() % MAX_J + 1;

	auto ref_mat =
		p_cuda.base_mut().alloc_ov_matrix<std::string, std::string>(
			I, J, meza::pool::pool_region::Reference);

	p_cuda.copy_region_to_host(meza::pool::pool_region::Reference);

	auto filter_mat =
		p_cuda.base_mut().alloc_ov_matrix<std::string, std::string>(
			I, J, meza::pool::pool_region::Filter);
	p_cuda.copy_region_to_host(meza::pool::pool_region::Filter);

	auto xor_mat =
		p_cuda.base_mut().alloc_ov_matrix<std::string, std::string>(
			I, J, meza::pool::pool_region::Xor);

	auto [ref_used, filter_used, xor_used] = p_cuda.base().used();
	p_cuda.sync();

	ASSERT_EQ(ref_used, I * J);
	ASSERT_EQ(filter_used, I * J);
	ASSERT_EQ(xor_used, I * J);

	ASSERT_FALSE(p_cuda.base().empty());
#endif
}

/* === alloc multi matrix tests === */

TEST(SplitPool, AllocCpuMultiMatrix)
{
	auto p = meza::pool::matrix_pool<u8>::create_from_megabytes();

	u32 MAX_I{10};
	u32 MAX_J{10};
	u32 MAX_MATRICES{10};

	u32 expected_used{};

	for (u32 ctr = 0; ctr < MAX_MATRICES; ctr++) {
		u32 I = std::rand() % MAX_I + 1;
		u32 J = std::rand() % MAX_J + 1;

		auto ref_mat = p.alloc_ov_matrix<std::string, std::string>(
			I, J, meza::pool::pool_region::Reference);

		auto filter_mat = p.alloc_ov_matrix<std::string, std::string>(
			I, J, meza::pool::pool_region::Filter);

		auto xor_mat = p.alloc_ov_matrix<std::string, std::string>(
			I, J, meza::pool::pool_region::Xor);

		auto [ref_used, filter_used, xor_used] = p.used();

		expected_used += I * J;

		ASSERT_EQ(ref_used, expected_used);
		ASSERT_EQ(filter_used, expected_used);
		ASSERT_EQ(xor_used, expected_used);
	}
}
