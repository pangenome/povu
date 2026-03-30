#ifndef IT_AT_MATRIX_NO_TANGLE_HPP
#define IT_AT_MATRIX_NO_TANGLE_HPP

#include <set> // for set

#include <liteseq/refs.h>      // for ref_walk, ref
#include <meza/pool/split.hpp> // for matrix_pool

#include "ita/traversals/at_matrix.hpp"	   // for rov_job_batch
#include "ita/traversals/depth_matrix.hpp" // for depth_matrix
#include "ita/variation/rov.hpp"	   // for RoV
#include "povu/common/core.hpp"		   // for pt

namespace ita::at_matrix::no_tangle
{
void from_no_tangle(const ir::RoV *rov,
		    const std::set<pt::u32> &to_call_ref_ids,
		    ita::depth_matrix::depth_matrix dm,
		    meza::pool::split::matrix_pool<qt::u8> &ov_pool,
		    rov_job_batch &batch);

} // namespace ita::at_matrix::no_tangle
#endif // IT_AT_MATRIX_NO_TANGLE_HPP
