#ifndef IT_UNTANGLE_HPP
#define IT_UNTANGLE_HPP

#include <string_view> // for string_view
#include <vector>      // for vector

#include "ita/genomics/allele.hpp" // for Exp, depth_matrix
#include "ita/variation/rov.hpp"   // for RoV

#include "povu/common/core.hpp"	     // for pt, idx_t, id_t, op_t
#include "povu/graph/bidirected.hpp" // for VG, bd

namespace ita::untangle
{
inline constexpr std::string_view MODULE = "povu::genomics::untangle";

using namespace ia;

std::vector<depth_matrix> untangle(const bd::VG &g,
				   const std::set<pt::u32> &to_call_ref_ids,
				   const depth_matrix &dm, const ir::RoV &rov);

race gen_race(const bd::VG &g, const std::vector<pt::id_t> &sorted_w,
	      pt::u32 h_idx);
} // namespace ita::untangle

#endif // IT_UNTANGLE_HPP
