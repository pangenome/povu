#include "ita/traversals/untangle.hpp"

#include <set>	   // for set, operator!=
#include <string>  // for basic_string, string
#include <utility> // for move
#include <vector>  // for vector

#include <quilt/graph_types.hpp> // for v_end_e, side_n_id_t, side_n_idx_t
#include <quilt/shim.hpp>	 // for contains
#include <quilt/types.hpp>	 // for qt

#include "povu/graph/bidirected.hpp" // for VG, bd
#include "povu/refs/refs.hpp"	     // for lq_strand_to_pv_or

#include "ita/align/align.hpp"	   // for align, aln_level_e
#include "ita/genomics/allele.hpp" // for Exp, itn_t
#include "ita/traversals/traversals.hpp" // for unroll_haps, itinerary, allele_traversal
					 //

namespace ita::traversals::untangle
{
namespace lq = liteseq;
using namespace ita::traversals::traversals;

ia::at_itn gen_at_itn(const bd::VG &g, const itinerary &steps_itn,
		      qt::u32 h_idx)
{
	std::vector<ptg::walk_t> allele_traversals;
	const lq::ref_walk *hw = g.get_ref_vec(h_idx)->walk; // the hap walk

	for (const allele_traversal &at : steps_itn) {
		ptg::walk_t w;
		for (qt::u32 step_idx : at) {
			ptg::or_e o =
				pr::lq_strand_to_pv_or(hw->strands[step_idx]);
			qt::u32 v_id = hw->v_ids[step_idx];
			w.emplace_back(ptg::id_or_t{v_id, o});
		}
		allele_traversals.emplace_back(std::move(w));
	}

	return ia::at_itn{std::move(allele_traversals)};
}

/**
 *
 */
void lineup(const std::string &edit_transcript, qt::u32 ref_h_idx,
	    qt::u32 ref_loop_count, qt::u32 alt_h_idx, qt::u32 alt_loop_count,
	    aln_chain &chain)
{
	qt::u32 i{}; // aln index
	qt::u32 j{}; // ref itn index
	qt::u32 k{}; // alt itn index

	const qt::u32 N = edit_transcript.size();
	const qt::u32 O = ref_loop_count;
	const qt::u32 P = alt_loop_count;

	for (; i < N && j < O && k < P; i++) {
		char c = edit_transcript[i];

		if (c == 'I')
			j++;
		else if (c == 'D')
			k++; // ignore that char (j) in ref
		else if (c == 'M' || c == 'X')
			chain.add({c, ref_h_idx, j++, alt_h_idx, k++});
		else
			throw std::runtime_error("Invalid edit script char");
	}
}

void do_align(const bd::VG &g, const std::set<qt::u32> &to_call_ref_ids,
	      aln_chain &chain)
{
	const qt::u32 I = g.get_hap_count();
	const std::vector<itinerary> &hap_itns = chain.hap_itns;

	// -----------
	// memoise fn
	//
	//
	// -----------
	std::map<qt::u32, ia::at_itn> hap2at_itn;
	auto get_at_itn = [&](qt::u32 h_idx) -> const ia::at_itn &
	{
		if (qs::contains(hap2at_itn, h_idx))
			return hap2at_itn[h_idx];

		const itinerary &itn = hap_itns[h_idx];
		hap2at_itn.emplace(h_idx, gen_at_itn(g, itn, h_idx));

		return hap2at_itn[h_idx];
	};

	// ------------------
	// run all alignments
	// ------------------
	const auto ALN_LEVEL = ita::align::aln_level_e::at;
	std::string et; // edit transcript
	for (qt::u32 ref_h_idx : to_call_ref_ids) {
		for (qt::u32 alt_h_idx{}; alt_h_idx < I; alt_h_idx++) {
			if (ref_h_idx == alt_h_idx)
				continue;

			const ia::at_itn &ref_itn = get_at_itn(ref_h_idx);
			const ia::at_itn &alt_itn = get_at_itn(alt_h_idx);

			et = ita::align::align(ref_itn, alt_itn, ALN_LEVEL);

			qt::u32 ref_loop_no = hap_itns[ref_h_idx].size();
			qt::u32 alt_loop_no = hap_itns[alt_h_idx].size();

			lineup(et, ref_h_idx, ref_loop_no, alt_h_idx,
			       alt_loop_no, chain);
		}

		// TODO: A: is this wise? theoretically
		if (!qs::contains(chain.all_chains, ref_h_idx)) {
			// when this is true entire alignment has no matches or
			// mismatches (only insertions and deletions), so we
			// skip it
			continue;
		}

		const std::map<qt::u32, std::vector<chain_link>> &loop2ats =
			chain.all_chains.at(ref_h_idx).loop2ats;

		// loop nos to remove from the chain
		std::set<qt::u32> to_remove;

		for (const auto &[loop_no, links] : loop2ats) {
			bool keep{false};
			for (auto &link : links)
				if (link.edit == 'X') {
					keep = true;
					break;
				}

			if (!keep)
				to_remove.insert(loop_no);
		}

		for (qt::u32 loop_no : to_remove)
			chain.all_chains.at(ref_h_idx).loop2ats.erase(loop_no);
	}
}

aln_chain untangle(const bd::VG &g, const std::set<qt::u32> &to_call_ref_ids,
		   const ir::RoV &rov)
{
	const std::vector<qt::u32> &sorted_w = rov.get_sorted_vertices();
	std::vector<itinerary> hap_itns = unroll_haps(g, sorted_w);
	aln_chain aln_chain{std::move(hap_itns)};
	do_align(g, to_call_ref_ids, aln_chain);

	return aln_chain;
}
} // namespace ita::traversals::untangle
