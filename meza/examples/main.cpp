#include <iostream>
#include <queue>
#include <string> // std::to_string
#include <vector>

#include <quilt/types.hpp>

#include "meza/pool/split.hpp"
#include "meza/shared/shared.hpp" // for layout

#include "./hap_xors.hpp"
#include "./shared.hpp"

using namespace shared;
using u32 = qt::u32;

void comp_conv(meza::pool::split::matrix_pool<qt::u8> &ov_pool, qt::u32 I,
	       qt::u32 total_J)
{
	// copy the pool to the device before launching the kernel
	ov_pool.copy_to_device();

	ov_pool.xor_on_device(total_J);

	// copy the pool back to the host after the kernel has finished
	ov_pool.copy_to_host_thirds(meza::pool::split::pool_region::Xor);
}

void copy_to_device(meza::pool::split::matrix_pool<qt::u8> &ov_pool)
{
#ifdef CONVO_USE_CUDA
	ov_pool.copy_to_device();
#endif
}

struct at_comparison {
	std::set<qt::up_t<u32>> matches;
	std::set<qt::up_t<u32>> mismatches;
};

/**
 *
 * compare allele traversals
 *
 * checks upper/super diagonals only
 */
at_comparison compare_ats(const std::vector<qt::u32> &ps, u32 mat_j, u32 I)
{
	std::set<qt::up_t<u32>> seen;
	std::set<qt::up_t<u32>> mismatches;
	std::set<qt::up_t<u32>> matches;

	u32 off_start = 0;

	const u32 K = I; // number of diagonals above the main diagonal

	std::cerr << __func__ << " mat_j " << mat_j << " I " << I << "\n";

	// r_idx is row index, c_idx is col index, k is diagonal
	for (u32 k{1}; k < K; k++) {
		for (u32 i{}; i < I; i++) { // row

			u32 j = i + k;
			if (i == j) // k == 0, main diagonal, skip
				continue;

			if (j >= I) // out of bounds, skip
				break;
			// diagonal & col
			// u32 c_idx = r_idx + k;

			if (seen.count({i, j}))
				continue;

			seen.insert({i, j});

			std::cerr << "comparing (" << i << ", " << j
				  << ") with k " << k << "\n";

			u32 off_end = (off_start + mat_j) - 1;

			std::cerr << "off idxs " << off_start << " " << off_end
				  << " vals " << ps[off_start] << " "
				  << ps[off_end] << "\n";

			u32 hap_a = i;
			u32 hap_b = j;

			bool start_is_match =
				(off_start == 0 ||
				 (off_start > 0 &&
				  ps[off_start] - ps[off_start - 1] == 0));

			if (start_is_match && ps[off_end] - ps[off_start] == 0)
				matches.insert({hap_a, hap_b});
			else
				mismatches.insert({hap_a, hap_b});

			off_start += mat_j;
		}
	}

	return {matches, mismatches};
}

// void comp_similaritites(meza::matrix_pool::split::matrix_pool<qt::u8>
// &ov_pool,
//			const std::vector<rov_mat_set> &mat_sets)
// {
//	auto to_u32 = [](qt::u8 v) -> u32
//	{
//		return static_cast<u32>(v);
//	};

//	for (u32 i{}; i < mat_sets.size(); i++) {

//		std::cerr << __func__ << " " << i << "\n";

//		const rov_mat_set &mat_set = mat_sets[i];
//		mat_set.dbg_print();

//		const auto &[_, filter_mat, ___, pool_j_offset] = mat_sets[i];
//		u32 mat_j = filter_mat.cols();

//		const u32 I = filter_mat.rows();
//		const u32 K = filter_mat.rows();

//		std::vector<qt::u8> res;
//		res.reserve(I * I / 2);

//		for (u32 k{1}; k < K; k++)
//			ov_pool.haps_xor(k, pool_j_offset, mat_j, I, res);

//		u32 N = res.size();
//		std::vector<qt::u32> ps(N, 0);

//		meza::matrix_pool::split::compare_rows(res, ps);

//		for (auto x : res)
//			std::cerr << to_u32(x) << ", ";
//		std::cerr << "\n";

//		std::cerr << "populated ps " << N << "\n";
//		for (auto x : ps)
//			std::cerr << to_u32(x) << ", ";
//		std::cerr << "\n";

//		const auto &[matches, mismatches] = compare_ats(ps, mat_j, I);

//		std::cerr << "Matches:\n";
//		for (const auto &[a, b] : matches)
//			std::cerr << "(" << a << ", " << b << ")\n";

//		std::cerr << "Mismatches:\n";
//		for (const auto &[a, b] : mismatches)
//			std::cerr << "(" << a << ", " << b << ")\n";

//		std::vector<qt::u8> res_sum;
//		res_sum.reserve(I * I / 2);
//		for (u32 k{1}; k < K; k++)
//			ov_pool.haps_sum(k, pool_j_offset, mat_j, I, res_sum);

//		// print the sums
//		std::cerr << "Sums:\n";
//		for (auto x : res_sum)
//			std::cerr << to_u32(x) << ", ";
//		std::cerr << "\n";
//	}
// }

int old_main()
{

	auto &ov_pool = meza::pool::split::matrix_pool<qt::u8>::init();

	std::queue<qt::u32> q;
	qt::u32 ctr = 0;
	qt::u32 max = 5;

	q.push(ctr);

	qt::u32 I = 7;
	qt::u32 J = 3;
	qt::u32 j_offset = 0;

	std::vector<meza::pool::split::rov_mat_set> mat_sets;

	while (!q.empty() && ctr < max) {
		std::cerr << "Iteration " << ctr << "\n";

		q.pop();

		if (!ov_pool.can_allocate(
			    I, J, meza::shared::layout::DenseRowMajor)) {
			comp_conv(ov_pool, I, j_offset);
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

		fill_random(filter_mat);
		fill_ref_row(filter_mat, ref_mat);

		mat_sets.push_back({ref_mat, filter_mat, xor_mat, j_offset});

		j_offset += J;
		q.push(ctr++);
	}

	qt::u32 expected_j_offset{};
	for (qt::u32 i{}; i < mat_sets.size(); i++) {
		const meza::pool::split::rov_mat_set &s = mat_sets[i];
		assert(s.j_offset == expected_j_offset);
		qt::u32 J = s.filter.cols();
		expected_j_offset += J;
	}

	std::cerr << __func__ << " offsets checked\n";

	std::cerr << ov_pool.used() << "\n";

	copy_to_device(ov_pool);
	// ov_pool.cuda_setup_haps_xor();

	if (!ov_pool.empty()) {
		// ov_pool.print_foo();
		ov_pool.copy_to_device();
		// comp_similaritites(ov_pool, mat_sets);
		//  ov_pool.sync_device();
		//  ov_pool.reset();
	}

	for (qt::u32 s{}; s < mat_sets.size(); s++) {
		std::cerr << "\n\n";
		mat_sets[s].dbg_print();
		std::cerr << "\n\n";
	}

	return 0;
}

int main()
{
#ifdef CONVO_DEBUG
	std::cerr << "Debug\n";
#elif CONVO_RELEASE
	std::cerr << "Release mode\n";
#endif

#ifdef CONVO_USE_CUDA

	std::cerr << "use\n";
#else
	std::cerr << "no use\n";
#endif

	auto &ov_pool = meza::pool::split::matrix_pool<qt::u8>::init();

	hap_xors::run(ov_pool);
	return 0;
}
