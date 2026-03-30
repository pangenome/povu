#ifndef IT_AT_MATRIX_TANGLED_HPP
#define IT_AT_MATRIX_TANGLED_HPP

#include <meza/pool/split.hpp> // for matrix_pool

#include "ita/traversals/at_matrix.hpp"
#include "ita/traversals/untangle.hpp"

namespace ita::at_matrix::tangled
{

void from_tangled(const bd::VG &g, const ir::RoV *rov,
		  const std::set<pt::u32> &to_call_ref_ids,
		  meza::pool::split::matrix_pool<qt::u8> &ov_pool,
		  const ita::traversals::untangle::aln_chain &ac,
		  rov_job_batch &batch);

} // namespace ita::at_matrix::tangled
#endif // IT_AT_MATRIX_TANGLED_HPP
