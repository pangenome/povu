#ifndef ITA_TRAVERSALS_HPP
#define ITA_TRAVERSALS_HPP

#include <vector> // for vector

#include <oza/graph/bidirected.hpp> // for VG, bd
#include <quilt/types.hpp>	    // for qt

namespace ita::traversals::traversals
{
using allele_traversal = std::vector<qt::u32>;
using itinerary = std::vector<allele_traversal>;

std::vector<itinerary> unroll_haps(const bd::VG &g,
				   const std::vector<qt::u32> &sorted_w);
} // namespace ita::traversals::traversals
#endif // ITA_TRAVERSALS_HPP
