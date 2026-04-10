#include <gtest/gtest.h>

#include <quilt/types.hpp>

#include <meza/pool/hap_comp.hpp>
#include <meza/pool/matrix_pool.hpp>
#include <meza/pool/split_pool_types.hpp>

#if MEZA_USE_CUDA
#include <meza/pool/hap_comp.cuh>
#include <meza/pool/matrix_pool_cuda.cuh>
#endif

using qt::u32;
using qt::u8;

using meza::pool::hap_comp::hap_comp_matrix;

TEST(HapComp, InitCpu)
{
	auto p = meza::pool::matrix_pool<u8>::create_from_megabytes();

	std::size_t capacity = 512 * 1024 * 1024;
	auto cmp_mat = hap_comp_matrix<u8>::create(capacity);

	u32 MAX_I{10};
	u32 MAX_J{10};

	u32 I = std::rand() % MAX_I + 1;
	u32 J = std::rand() % MAX_J + 1;

	auto ref_mat = p.alloc_ov_matrix<std::string, std::string>(
		I, J, meza::pool::pool_region::Reference);

	auto filter_mat = p.alloc_ov_matrix<std::string, std::string>(
		I, J, meza::pool::pool_region::Filter);

	cmp_mat.set_filter(&filter_mat, 0);

	auto xor_mat = p.alloc_ov_matrix<std::string, std::string>(
		I, J, meza::pool::pool_region::Xor);

	auto [ref_used, filter_used, xor_used] = p.used();

	ASSERT_EQ(ref_used, I * J);
	ASSERT_EQ(filter_used, I * J);
	ASSERT_EQ(xor_used, I * J);

	ASSERT_FALSE(p.empty());
}

TEST(HapComp, InitCuda)
{
#if MEZA_USE_CUDA
	auto p = meza::pool::matrix_pool<u8>::create_from_megabytes();

	std::size_t capacity = 512 * 1024 * 1024;
	auto cmp_mat = hap_comp_matrix<u8>::create(capacity);
	auto cmp_mat_cuda =
		meza::pool::hap_comp::hap_comp_matrix_cuda<u8>{cmp_mat};

	u32 MAX_I{10};
	u32 MAX_J{10};

	u32 I = std::rand() % MAX_I + 1;
	u32 J = std::rand() % MAX_J + 1;

	auto ref_mat = p.alloc_ov_matrix<std::string, std::string>(
		I, J, meza::pool::pool_region::Reference);

	auto filter_mat = p.alloc_ov_matrix<std::string, std::string>(
		I, J, meza::pool::pool_region::Filter);

	cmp_mat_cuda.base_mut().set_filter(&filter_mat, 0);

	auto xor_mat = p.alloc_ov_matrix<std::string, std::string>(
		I, J, meza::pool::pool_region::Xor);

	auto [ref_used, filter_used, xor_used] = p.used();

	ASSERT_EQ(ref_used, I * J);
	ASSERT_EQ(filter_used, I * J);
	ASSERT_EQ(xor_used, I * J);

	ASSERT_FALSE(p.empty());
#endif
}
