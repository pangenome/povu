#include "ita/traversals/traversals.hpp"

#include <algorithm>
#include <utility> // for move
#include <vector>  // for vector

#include <quilt/types.hpp> // for qt

#include "povu/graph/bidirected.hpp" // for VG, bd

namespace ita::traversals::traversals
{
/**
 * assumption: no duplicate steps can occur
 */
std::vector<qt::u32> lineup(const bd::VG &g,
			    const std::vector<qt::u32> &sorted_w, qt::u32 h_idx)
{
	std::vector<qt::u32> unrolled;
	for (qt::u32 j{}; j < sorted_w.size(); j++) {
		qt::u32 v_idx = g.v_id_to_idx(sorted_w[j]);
		for (qt::u32 step_idx : g.get_vertex_ref_idxs(v_idx, h_idx))
			unrolled.push_back(step_idx);
	}

	std::sort(unrolled.begin(), unrolled.end());

	return unrolled;
}

itinerary cluster(const std::vector<qt::u32> &unrolled)
{
	itinerary itn;
	if (unrolled.empty())
		return itn;

	// create a buffer and init with the first element from unrolled
	allele_traversal buf = {unrolled[0]};

	for (std::size_t i{1}; i < unrolled.size(); ++i) {
		if (unrolled[i - 1] + 1 == unrolled[i]) {
			buf.push_back(unrolled[i]);
		}
		else {
			itn.push_back(std::move(buf));
			buf = allele_traversal{};
			buf.push_back(unrolled[i]);
		}
	}

	itn.push_back(std::move(buf)); // flush remaining elements from buffer
	return itn;
}

std::vector<itinerary> unroll_haps(const bd::VG &g,
				   const std::vector<qt::u32> &sorted_w)
{
	std::vector<itinerary> hap_itns;

	qt::u32 I = g.get_hap_count();
	for (qt::u32 h_idx{}; h_idx < I; h_idx++) {
		std::vector<qt::u32> unrolled = lineup(g, sorted_w, h_idx);
		itinerary itn = cluster(unrolled);
		hap_itns.emplace_back(std::move(itn));
	}

	return hap_itns;
}
} // namespace ita::traversals::traversals
