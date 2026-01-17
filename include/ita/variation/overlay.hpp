#ifndef IT_OVERLAY_HPP
#define IT_OVERLAY_HPP

#include <set>	  // for set
#include <vector> // for vector

#include "ita/genomics/allele.hpp" // for Exp
#include "ita/variation/rov.hpp"   // for RoV
#include "ita/variation/sne.hpp"   // for pin_cushion

namespace ita::overlay
{
inline constexpr std::string_view MODULE = "povu::overlay";

namespace lq = liteseq;
namespace pgt = povu::types::graph;
namespace pvst = povu::pvst;

constexpr ir::var_type_e ins = ir::var_type_e::ins;
constexpr ir::var_type_e del = ir::var_type_e::del;
constexpr ir::var_type_e sub = ir::var_type_e::sub;
constexpr pgt::or_e fo = pgt::or_e::forward;
constexpr pgt::or_e ro = pgt::or_e::reverse;

std::vector<ia::trek> overlay_generic(const bd::VG &g, ir::RoV &rov,
				      const std::set<pt::u32> &to_call_ref_ids,
				      ise::pin_cushion &pcushion);

} // namespace ita::overlay

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace po = ita::overlay;

#endif // IT_OVERLAY_HPP
