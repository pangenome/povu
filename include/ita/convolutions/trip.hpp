#ifndef ITA_TRIP_HPP
#define ITA_TRIP_HPP

#include <convo/pool.hpp> // for matrix_pool
#include <liteseq/refs.h> // for ref_walk, ref

#include "ita/convolutions/at_matrix.hpp" // for matrix_pool, rov_matrix_set
#include "quilt/types.hpp"		  // for qt::u32

namespace ita::trip
{

struct at_comparison {
	std::set<qt::up_t<qt::u32>> matches;
	std::set<qt::up_t<qt::u32>> mismatches;
};

at_comparison
comp_at_comparison(meza::matrix_pool::matrix_pool<qt::u8> &ov_pool,
		   const ita::at_matrix::mat3 &mat_set);
}; // namespace ita::trip

#endif // ITA_TRIP_HPP
