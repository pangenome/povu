#ifndef ITA_CONVOLUTIONS_HPP
#define ITA_CONVOLUTIONS_HPP

#include <set>
#include <vector>

#include <convo/pool.hpp> // for matrix_pool

#include "ita/convolutions/at_matrix.hpp" // rov_matrix_set
#include "ita/genomics/allele.hpp"	  // for trek
#include "ita/variation/rov.hpp"	  // for RoV
#include "povu/common/core.hpp"
#include "povu/graph/bidirected.hpp" // for VG

namespace ita::convolutions
{

void populate_trips(const bd::VG &g, const ir::RoV *rov,
		    const std::set<pt::u32> &to_call_ref_ids,
		    meza::matrix_pool::joint_pool<qt::u32> &dm_pool,
		    meza::matrix_pool::matrix_pool<qt::u8> &ov_pool,
		    ita::at_matrix::rov_job_batch &batch,
		    std::vector<ia::trek> &treks, bool drain);

// void comp_expedition(const bd::VG &g, const ir::RoV *rov, bool last,
//		     const std::set<pt::u32> &to_call_ref_ids,
//		     meza::matrix_pool::matrix_pool<qt::u8> &ov_pool,
//		     ita::at_matrix::all_mat_sets &mat_sets,
//		     std::vector<ia::trek> &treks);
} // namespace ita::convolutions

#endif // ITA_CONVOLUTIONS_HPP
