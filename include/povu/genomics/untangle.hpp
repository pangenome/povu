#ifndef POVU_UNTANGLE_HPP
#define POVU_UNTANGLE_HPP

#include <string_view> // for string_view
#include <vector>

#include "povu/common/core.hpp"
#include "povu/genomics/allele.hpp" // for Exp

namespace povu::genomics::untangle
{
inline constexpr std::string_view MODULE = "povu::genomics::untangle";
namespace pga = povu::genomics::allele;
using namespace povu::genomics::allele;

std::vector<depth_matrix> untangle(const bd::VG &g,
				   const std::set<pt::u32> &to_call_ref_ids,
				   const depth_matrix &dm, const pvr::RoV &rov);

race gen_race(const bd::VG &g, const std::vector<pt::id_t> &sorted_w,
	      pt::u32 h_idx);
} // namespace povu::genomics::untangle

#endif // POVU_UNTANGLE_HPP
