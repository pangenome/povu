#ifndef IT_AT_MATRIX_TANGLED_HPP
#define IT_AT_MATRIX_TANGLED_HPP

#include "ita/convolutions/at_matrix.hpp"
#include "ita/traversals/untangle.hpp"

namespace ita::at_matrix::tangled
{
ita::at_matrix::rov_matrix_pool
init_tangled_depth_matrices(const bd::VG &g, const ir::RoV &rov,
			    const std::set<pt::u32> &to_call_ref_ids,
			    const ita::traversals::untangle::aln_chain &c);
} // namespace ita::at_matrix::tangled
#endif // IT_AT_MATRIX_TANGLED_HPP
