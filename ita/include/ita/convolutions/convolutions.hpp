#ifndef ITA_CONVOLUTIONS_HPP
#define ITA_CONVOLUTIONS_HPP

#include <set>
#include <string_view>
#include <vector>

#include <meza/pool/joint.hpp> // for joint_pool
#include <meza/pool/split.hpp> // for matrix_pool

#include "ita/genomics/allele.hpp"	// for trek
#include "ita/traversals/at_matrix.hpp" // rov_matrix_set
#include "ita/variation/rov.hpp"	// for RoV
#include "povu/common/core.hpp"
#include "povu/graph/bidirected.hpp" // for VG

namespace ita::convolutions
{

constexpr std::string_view MODULE = "ita::convolutions";

void populate_trips(const bd::VG &g, const ir::RoV *rov,
		    const std::set<pt::u32> &to_call_ref_ids,
		    meza::pool::joint::joint_pool<qt::u32> &dm_pool,
		    meza::pool::split::matrix_pool<qt::u8> &ov_pool,
		    ita::at_matrix::rov_job_batch &batch,
		    std::vector<ia::trek> &treks, bool drain);
} // namespace ita::convolutions

#endif // ITA_CONVOLUTIONS_HPP
