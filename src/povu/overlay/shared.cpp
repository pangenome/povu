#include "povu/overlay/overlay.hpp"

#include <liteseq/refs.h> // for ref_walk, ref
#include <map>		  // for map
#include <utility>

#include "povu/common/core.hpp"
#include "povu/genomics/allele.hpp"
#include "povu/overlay/shared.hpp"

namespace povu::overlay::shared
{
void update_exp(pt::u32 r_idx, pt::u32 w_idx, pga::allele_slice_t &&at,
		pga::Exp &e)
{
	std::map<pt::id_t, pga::itn_t> &ref_map = e.get_ref_itns_mut();
	pga::itn_t &itn = ref_map[r_idx];

	std::map<pt::idx_t, std::set<pt::idx_t>> &walk_to_refs =
		e.get_walk_to_ref_idxs_mut();

	itn.append_at_sorted(std::move(at));
	walk_to_refs[w_idx].insert(r_idx);
}
} // namespace povu::overlay::shared
