#ifndef IT_AT_MATRIX_TANGLED_HPP
#define IT_AT_MATRIX_TANGLED_HPP

#include <meza/pool/pool.hpp> // for pool

#include "ita/traversals/at_matrix.hpp"
#include "ita/traversals/untangle.hpp"

namespace ita::at_matrix::tangled
{
using pool_t = meza::pool::pool<qt::u8, qt::u32>;

void from_tangled(const bd::VG &g, const ir::RoV *rov,
		  const std::set<pt::u32> &to_call_ref_ids, pool_t &p,
		  const ita::traversals::untangle::aln_chain &ac,
		  rov_job_batch &batch);

} // namespace ita::at_matrix::tangled
#endif // IT_AT_MATRIX_TANGLED_HPP
