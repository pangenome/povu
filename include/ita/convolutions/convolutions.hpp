#ifndef ITA_CONVOLUTIONS_HPP
#define ITA_CONVOLUTIONS_HPP

#include <optional>
#include <set>
#include <vector>

#include "ita/genomics/allele.hpp" // for trek
#include "ita/variation/rov.hpp"   // for RoV
#include "povu/common/core.hpp"
#include "povu/graph/bidirected.hpp" // for VG

namespace ita::convolutions
{

std::optional<ia::trek>
comp_expedition(const bd::VG &g, ir::RoV &rov,
		const std::set<pt::u32> &to_call_ref_ids);
} // namespace ita::convolutions

#endif // ITA_CONVOLUTIONS_HPP
