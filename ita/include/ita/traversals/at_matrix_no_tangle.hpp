#ifndef IT_AT_MATRIX_NO_TANGLE_HPP
#define IT_AT_MATRIX_NO_TANGLE_HPP

#include <set> // for set

#include <liteseq/refs.h>     // for ref_walk, ref
#include <meza/pool/pool.hpp> // for matrix_pool
#include <quilt/types.hpp>    // for qt

#include "ita/traversals/at_matrix.hpp"	   // for rov_job_batch
#include "ita/traversals/depth_matrix.hpp" // for depth_matrix
#include "ita/variation/rov.hpp"	   // for RoV

namespace ita::at_matrix::no_tangle
{
using pool_t = meza::pool::pool<qt::u8, qt::u32>;

void from_no_tangle(const ir::RoV *rov,
		    const std::set<qt::u32> &to_call_ref_ids,
		    ita::depth_matrix::depth_matrix dm, pool_t &p,
		    rov_job_batch &batch);

} // namespace ita::at_matrix::no_tangle
#endif // IT_AT_MATRIX_NO_TANGLE_HPP
