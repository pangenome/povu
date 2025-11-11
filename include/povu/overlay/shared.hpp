#ifndef POVU_OVERLAY_SHARED_HPP
#define POVU_OVERLAY_SHARED_HPP

// #include "povu/common/constants.hpp"
#include "povu/common/core.hpp"	    // for pt, idx_t, id_t, op_t
#include "povu/genomics/allele.hpp" // for Exp

namespace povu::overlay::shared
{
void update_exp(pt::u32 r_idx, pt::u32 w_idx, pga::allele_slice_t &&at,
		pga::Exp &e);
}

#endif // POVU_OVERLAY_SHARED_HPP
