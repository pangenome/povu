#ifndef ITA_TRAVERSALS_HPP
#define ITA_TRAVERSALS_HPP

#include <vector> // for vector

// #include "ita/convolutions/at_matrix.hpp" // for matrix_pool
#include "povu/common/core.hpp"	     // for pt
#include "povu/graph/bidirected.hpp" // for VG, bd

namespace ita::traversals::traversals
{
using allele_traversal = std::vector<pt::u32>;
using itinerary = std::vector<allele_traversal>;

std::vector<itinerary> unroll_haps(const bd::VG &g,
				   const std::vector<pt::u32> &sorted_w);
} // namespace ita::traversals::traversals
#endif // ITA_TRAVERSALS_HPP
