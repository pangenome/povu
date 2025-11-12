#include "povu/overlay/overlay.hpp"

#include <liteseq/refs.h> // for ref_walk, ref
#include <utility>
#include <vector> // for vector

#include "povu/common/core.hpp"
#include "povu/genomics/allele.hpp"
#include "povu/overlay/generic.hpp"
#include "povu/overlay/sne.hpp"
#include "povu/variation/rov.hpp"

namespace povu::overlay
{
namespace lq = liteseq;

std::pair<std::vector<pga::Exp>, std::vector<pos::pin_cushion>>
comp_itineraries3(const bd::VG &g, const pvr::RoV &rov,
		  const std::set<pt::id_t> &to_call_ref_ids)
{
	if (rov.can_be_non_planar())
		generic::overlay_generic(g, rov);

	auto [e, opt_pc] = overlay_tiny(g, rov, to_call_ref_ids);
	if (opt_pc.has_value())
		return {{e}, {std::move(opt_pc.value())}};
	else
		return {{e}, {}};
}
} // namespace povu::overlay
