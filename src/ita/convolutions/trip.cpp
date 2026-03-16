#include "ita/convolutions/trip.hpp"

// #include <iostream> // for std::cerr
#include <set>	  // for std::set
#include <string> // for std::string
#include <vector> // for std::vector

#include <convo/pool.hpp> // for matrix_pool
#include <liteseq/refs.h> // for ref_walk, ref

#include "ita/convolutions/at_matrix.hpp" // for matrix_pool, rov_matrix_set
#include "quilt/types.hpp"

namespace ita::trip
{

using ov_mat_t = meza::matrix_view::ov_matrix<qt::u8, std::string, std::string>;
using u32 = qt::u32;

u32 u8_to_u32(qt::u8 v)
{
	return static_cast<u32>(v);
};

/**
 *
 * compare allele traversals
 *
 * checks upper/super diagonals only
 */
at_comparison compare_ats(const ov_mat_t &filter_mat,
			  const std::vector<qt::u32> &ps, u32 mat_j, u32 I)
{
	std::set<qt::up_t<u32>> seen;

	std::set<qt::up_t<u32>> mismatches;
	std::set<qt::up_t<u32>> matches;

	u32 off_start = 0;

	const u32 K = I; // number of diagonals above the main diagonal

	auto is_blank = [&](u32 i, u32 j) -> bool
	{
		return filter_mat.base().is_row_blank(i) ||
		       filter_mat.base().is_row_blank(j);
	};

	// std::cerr << __func__ << " mat_j " << mat_j << " I " << I << "\n";

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

			if (is_blank(i, j)) {
				off_start += mat_j;
				continue;
			}

			if (seen.count({i, j}))
				continue;

			seen.insert({i, j});

			// std::cerr << "comparing (" << i << ", " << j
			//	  << ") with k " << k << "\n";

			u32 off_end = (off_start + mat_j) - 1;

			// std::cerr << "off idxs " << off_start << " " <<
			// off_end
			//	  << " vals " << ps[off_start] << " "
			//	  << ps[off_end] << "\n";

			// for (qt::u32 i{off_start}; i < off_end; i++)
			//	std::cerr << u8_to_u32(ps[i]) << ", ";
			// std::cerr << "\n";

			u32 hap_a = i;
			u32 hap_b = j;

			// check the one before offs start because it
			// has the sum
			off_start > 0 ? off_start-- : off_start;

			// bool start_is_match =
			//	(off_start == 0 ||
			//	 (off_start > 0 &&
			//	  ps[off_start] - ps[off_start - 1] == 0));

			if (ps[off_end] - ps[off_start] == 0)
				matches.insert({hap_a, hap_b});
			else
				mismatches.insert({hap_a, hap_b});

			off_start += mat_j;
		}
	}

	std::cerr << "matches:\n";
	for (auto [a, b] : matches) {
		std::cerr << "(" << a << ", " << b << ")\n";
	}

	std::cerr << "mismatches:\n";
	for (auto [a, b] : mismatches) {
		std::cerr << "(" << a << ", " << b << ")\n";
	}

	return {matches, mismatches};
}

at_comparison
comp_at_comparison(meza::matrix_pool::matrix_pool<qt::u8> &ov_pool,
		   const ita::at_matrix::mat3 &mat_set)
{

	// mat_set.dbg_print();
	const ov_mat_t &filter_mat = mat_set.filter;
	qt::u32 pool_j_offset = mat_set.j_offset;

	// const auto &[_, filter_mat, __, pool_j_offset, ___] = mat_set;
	u32 mat_j = filter_mat.cols();
	const u32 I = filter_mat.rows();
	const u32 K = I;

	std::vector<qt::u8> res;
	res.reserve(I * I / 2);

	for (u32 k{1}; k < K; k++)
		ov_pool.haps_xor(k, pool_j_offset, mat_j, I, res);

	u32 N = res.size();
	std::vector<qt::u32> ps(N, 0);

	meza::matrix_pool::compare_rows(res, ps);

	// print_rows(filter_mat, res, mat_j, I);
	// std::cerr << "populated ps " << N << "\n";
	// for (auto x : ps)
	//	std::cerr << to_u32(x) << ", ";
	// std::cerr << "\n";

	at_comparison comp = compare_ats(filter_mat, ps, mat_j, I);

	// std::cerr << "Matches:\n";
	// for (const auto &[a, b] : matches)
	//	std::cerr << "(" << a << ", " << b << ")\n";

	// std::cerr << "Mismatches:\n";
	// for (const auto &[a, b] : mismatches)
	//	std::cerr << "(" << a << ", " << b << ")\n";

	std::vector<qt::u8> res_sum;
	res_sum.reserve(I * I / 2);
	for (u32 k{1}; k < K; k++)
		ov_pool.haps_sum(k, pool_j_offset, mat_j, I, res_sum);

	// print the sums
	// std::cerr << "Sums:\n";
	// for (auto x : res_sum)
	//	std::cerr << to_u32(x) << ", ";
	// std::cerr << "\n";

	return comp;
}

}; // namespace ita::trip
