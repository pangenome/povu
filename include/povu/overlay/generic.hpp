#ifndef POVU_OVERLAY_GENERIC_HPP
#define POVU_OVERLAY_GENERIC_HPP

#include <set>	   // for set
#include <utility> // for pair
#include <vector>  // for vector

// #include "povu/common/constants.hpp"
#include "povu/common/core.hpp"	    // for pt, idx_t, id_t, op_t
#include "povu/genomics/allele.hpp" // for Exp
#include "povu/variation/rov.hpp"   // for RoV

namespace povu::overlay::generic
{
std::pair<std::vector<pga::Exp>, std::vector<pga::sub_inv>>
overlay_generic(const bd::VG &g, const pvr::RoV &rov);
};
#endif // POVU_OVERLAY_GENERIC_HPP
