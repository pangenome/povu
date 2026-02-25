#ifndef ITA_AT_MATRIX_HPP
#define ITA_AT_MATRIX_HPP

#include <map> // for map
#include <set> // for set

#include <liteseq/refs.h> // for ref_walk, ref

#include "ita/genomics/allele.hpp"	 // for hap_slice
#include "ita/traversals/traversals.hpp" // for itinerary
#include "ita/variation/rov.hpp"	 // for RoV
#include "meza/matrix/matrix.hpp"	 // for matrix2d
#include "povu/common/core.hpp"		 // for pt
#include "povu/graph/bidirected.hpp"	 // for VG

namespace ita::at_matrix
{
struct matrix_pool {
	std::map<pt::u32, meza::matrix::ref_matrix> ref_matrices;
	meza::matrix::depth_matrix filter_matrix;
	meza::matrix::depth_matrix result_matrix;
	std::vector<std::optional<ia::hap_slice>> hap_slices;

	// only used for tangled pools, otherwise 0 by default
	pt::u32 ref_loop_no = 0;

	bool tangled;
	pt::u32 I;
	pt::u32 J;
};

struct rov_matrix_pool {
	const ir::RoV &rov;
	std::vector<matrix_pool> pools;
	bool tangled = false; // if not tangled pools should have size 1
	std::vector<ita::traversals::traversals::itinerary> hap_itns;

	// --------------
	// constructor(s)
	// --------------

	rov_matrix_pool(const ir::RoV &rov) : rov(rov)
	{}

	rov_matrix_pool(
		const ir::RoV &rov, bool t,
		const std::vector<ita::traversals::traversals::itinerary>
			&hap_itns)
	    : rov(rov), tangled(t), hap_itns(hap_itns)
	{}

	// ---------
	// getter(s)
	// ---------

	[[nodiscard]]
	bool empty() const
	{
		return pools.empty();
	}
};

rov_matrix_pool init_depth_matrices(const bd::VG &g, ir::RoV &rov,
				    const std::set<pt::u32> &to_call_ref_ids);
} // namespace ita::at_matrix
#endif // ITA_AT_MATRIX_HPP
