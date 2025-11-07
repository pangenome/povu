#ifndef POVU_OVERLAY_HPP
#define POVU_OVERLAY_HPP

#include <set>	   // for set
#include <utility> // for pair
#include <vector>  // for vector

// #include "povu/common/constants.hpp"
#include "povu/common/core.hpp"	    // for pt, idx_t, id_t, op_t
#include "povu/genomics/allele.hpp" // for Exp
#include "povu/overlay/sne.hpp"	    // for pin_cushion
#include "povu/variation/rov.hpp"   // for RoV

namespace povu::overlay
{
inline constexpr std::string_view MODULE = "povu::overlay";

namespace lq = liteseq;
namespace pgt = povu::types::graph;
namespace pvst = povu::pvst;

constexpr pvr::var_type_e ins = pvr::var_type_e::ins;
constexpr pvr::var_type_e del = pvr::var_type_e::del;
constexpr pvr::var_type_e sub = pvr::var_type_e::sub;
constexpr pgt::or_e fo = pgt::or_e::forward;
constexpr pgt::or_e ro = pgt::or_e::reverse;

std::pair<std::vector<pga::Exp>, std::vector<pos::pin_cushion>>
comp_itineraries3(const bd::VG &g, const pvr::RoV &rov,
		  const std::set<pt::id_t> &to_call_ref_ids);

std::pair<pga::Exp, std::optional<pos::pin_cushion>>
overlay_tiny(const bd::VG &g, const pvr::RoV &rov,
	     const std::set<pt::id_t> &to_call_ref_ids);

} // namespace povu::overlay

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace po = povu::overlay;

#endif // POVU_OVERLAY_HPP
