#ifndef ITA_TR_UN_HPP
#define ITA_TR_UN_HPP

#include <map>
// #include <string>
#include <vector>

#include "ita/traversals/traversals.hpp" // for itinerary
#include "ita/variation/rov.hpp"	 // for RoV
#include "povu/common/core.hpp"		 // for pt

namespace ita::traversals::untangle
{
struct chain_link {
	char edit; // 'M', 'I', 'D'
	pt::u32 ref_h_idx;
	pt::u32 ref_loop_no;
	pt::u32 alt_h_idx;
	pt::u32 alt_loop_no;
};

struct chain {
	// loop no to aligned_ats
	std::map<pt::u32, std::vector<chain_link>> loop2ats;

	// ---------
	// getter(s)
	// ---------

	void add(chain_link &&link)
	{
		this->loop2ats[link.ref_loop_no].emplace_back(link);
	}
};

struct aln_chain {
	// ref_h_idx to chain
	std::map<pt::u32, chain> all_chains;
	std::vector<ita::traversals::traversals::itinerary> hap_itns;
	std::map<pt::u32, pt::u32> ref2max_loop_no;

	// --------------
	// constructor(s)
	// --------------

	aln_chain(
		std::vector<ita::traversals::traversals::itinerary> &&hap_itns)
	    : hap_itns(std::move(hap_itns))
	{}

	// ---------
	// getter(s)
	// ---------

	[[nodiscard]]
	pt::u32 get_max_loop_no(pt::u32 ref_h_idx) const
	{
		return this->ref2max_loop_no.at(ref_h_idx);
	}

	[[nodiscard]]
	std::optional<chain_link> get_by(pt::u32 ref_h_idx, pt::u32 alt_h_idx,
					 pt::u32 loop_no) const
	{
		const std::vector<chain_link> &cl =
			this->all_chains.at(ref_h_idx).loop2ats.at(loop_no);

		for (auto &link : cl)
			if (link.alt_h_idx == alt_h_idx)
				return link;

		return std::nullopt;
	}

	// ---------
	// setter(s)
	// ---------

	void add(const chain_link &link)
	{
		pt::u32 ref_h_idx = link.ref_h_idx;
		pt::u32 loop_no = link.ref_loop_no;

		// initialises to 0 if ref_h_idx not present
		ref2max_loop_no[ref_h_idx] =
			std::max(ref2max_loop_no[ref_h_idx], loop_no);

		chain_link cpy = link;
		this->all_chains[link.ref_h_idx].add(std::move(cpy));
	}
};

aln_chain untangle(const bd::VG &g, const std::set<pt::u32> &to_call_ref_ids,
		   const ir::RoV &rov);

} // namespace ita::traversals::untangle
#endif
