#ifndef ITA_TRIP_HPP
#define ITA_TRIP_HPP

#include <liteseq/refs.h>	  // for ref_walk, ref
#include <meza/pool/hap_comp.hpp> // for haps_comp_set
// #include <meza/pool/split.hpp>	  // for matrix_pool

#include <quilt/types.hpp> // for qt

#include "ita/traversals/at_matrix.hpp"	 // for matrix_pool, rov_matrix_set
#include "ita/traversals/traversals.hpp" // for itinerary
#include "quilt/types.hpp"		 // for qt::u32

namespace ita::trip
{

std::optional<ia::trek>
gen_trip(const bd::VG &g, const ir::RoV *rov, bool is_tangled,
	 qt::u32 ref_h_idx, const ita::at_matrix::hap2loop &h2l,
	 const std::vector<ita::traversals::traversals::itinerary> &hap_itns,
	 const ita::at_matrix::mat3 &mat_set,
	 const std::vector<qt::u32> &sorted_vertices,
	 const meza::pool::hap_comp::haps_comp_set &hap_cmp);

}; // namespace ita::trip

#endif // ITA_TRIP_HPP
