#include "povu/genomics/untangle.hpp"

#include <set>	   // for set, operator!=
#include <string>  // for basic_string, string
#include <utility> // for move
#include <vector>  // for vector

#include "fmt/core.h"		    // for format
#include "povu/align/align.hpp"	    // for align, aln_level_e
#include "povu/common/compat.hpp"   // for format, pv_cmp
#include "povu/common/core.hpp"	    // for pt, id_t, up_t, operator<
#include "povu/genomics/allele.hpp" // for Exp, itn_t

namespace povu::genomics::untangle
{
namespace pa = povu::align;

inline std::vector<pt::up_t<pt::id_t>> compute_pairs(const pga::Exp &rt)
{
	std::set<pt::id_t> ref_ids = rt.get_ref_ids();

	std::set<pt::up_t<pt::id_t>> done;
	std::vector<pt::up_t<pt::id_t>> aln_pairs;

	// all vs all compare
	for (pt::id_t ref_id1 : ref_ids) {
		for (pt::id_t ref_id2 : ref_ids) {

			if (ref_id1 == ref_id2) {
				continue;
			}

			pt::up_t<pt::id_t> p{ref_id1, ref_id2};

			if (done.find(p) != done.end()) {
				continue;
			}

			done.insert(p);
			aln_pairs.push_back(p);
		}
	}

	return aln_pairs;
}

// change I to D and D to I
std::string invert_aln(const std::string &aln)
{
	std::string inverted_aln;
	inverted_aln.reserve(aln.size());

	for (char c : aln) {
		switch (c) {
		case 'I':
			inverted_aln += 'D';
			break;
		case 'D':
			inverted_aln += 'I';
			break;
		default:
			inverted_aln += c; // keep other characters unchanged
		}
	}

	return inverted_aln;
}

/**
 * formerly untangle_flb
 * align the traversals of two refs
 */
void untangle_ref_walks(pga::Exp &rt)
{
	std::string et;

	for (pt::id_t ref_id1 : rt.get_ref_ids()) {
		for (pt::id_t ref_id2 : rt.get_ref_ids()) {
			if (ref_id1 == ref_id2) {
				continue;
			}

			if (rt.has_aln(ref_id2, ref_id1)) {
				// invert aln
				et = invert_aln(rt.get_aln(ref_id2, ref_id1));
			}
			else {
				const pga::itn_t &itn1 = rt.get_itn(ref_id1);
				const pga::itn_t &itn2 = rt.get_itn(ref_id2);

				et = pa::align(itn1, itn2, pa::aln_level_e::at);
			}

			rt.add_aln(ref_id1, ref_id2, std::move(et));
		}
	}

	std::vector<pt::up_t<pt::id_t>> aln_pairs = compute_pairs(rt);

	// for (auto [ref_id1, ref_id2] : aln_pairs) {

	//   const pga::Itn &itn1 = rt.get_itn(ref_id1);
	//   const pga::Itn &itn2 = rt.get_itn(ref_id2);

	//   std::string et = pa::align(itn1, itn2, pvt::aln_level_e::at);

	//   rt.add_aln(ref_id1, ref_id2, std::move(et));
	// }
}

} // namespace povu::genomics::untangle
