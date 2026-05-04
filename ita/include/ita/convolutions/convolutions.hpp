#ifndef ITA_CONVOLUTIONS_HPP
#define ITA_CONVOLUTIONS_HPP

#include <set>
#include <string_view>
#include <vector>

#include <meza/pool/pool.hpp>	    // for matrix_pool
#include <oza/graph/bidirected.hpp> // for VG
#include <quilt/types.hpp>	    // for qt

#include "ita/genomics/allele.hpp"	// for trek
#include "ita/traversals/at_matrix.hpp" // rov_matrix_set
#include "ita/variation/rov.hpp"	// for RoV

namespace ita::convolutions
{
constexpr std::string_view MODULE = "ita::convolutions";

using pool_t = meza::pool::pool<qt::u8, qt::u32>;

void populate_trips(const bd::VG &g, const ir::RoV *rov,
		    const std::set<qt::u32> &to_call_ref_ids, pool_t &p,
		    ita::at_matrix::rov_job_batch &batch,
		    std::vector<ia::trek> &treks, bool drain);
} // namespace ita::convolutions

#endif // ITA_CONVOLUTIONS_HPP
