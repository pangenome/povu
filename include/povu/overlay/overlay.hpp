#ifndef POVU_OVERLAY_HPP
#define POVU_OVERLAY_HPP

#include <set>	   // for set
#include <utility> // for pair
#include <vector>  // for vector

#include "povu/common/core.hpp"	    // for pt, idx_t, id_t, op_t
#include "povu/genomics/allele.hpp" // for Exp
#include "povu/variation/rov.hpp"   // for RoV

namespace povu::overlay
{
inline constexpr std::string_view MODULE = "povu::overlay";

namespace lq = liteseq;

namespace pgt = povu::types::graph;
namespace pvst = povu::pvst;

struct sub_inv {
	pt::u32 walk_idx;
	std::vector<pt::u32> fwd_refs;
	std::vector<pt::u32> rev_refs;
};

std::pair<std::vector<pga::Exp>, std::vector<sub_inv>>
comp_itineraries3(const bd::VG &g, const pvr::RoV &rov,
		  const std::set<pt::id_t> &to_call_ref_ids);

std::pair<std::vector<pga::Exp>, std::vector<sub_inv>>
overlay_tiny(const bd::VG &g, const pvr::RoV &rov,
	     const std::set<pt::id_t> &to_call_ref_ids);

} // namespace povu::overlay

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace po = povu::overlay;

#endif // POVU_OVERLAY_HPP
