#include "zien/common/common.hpp"

#include <set>
#include <string>

#include <liteseq/refs.h> // for ref_walk, ref

#include "mto/from_vcf.hpp"	     // for VCFile
#include "povu/common/core.hpp"	     // for pt
#include "povu/graph/bidirected.hpp" // for bidirected

namespace zien::common
{

char to_char(liteseq::strand s)
{
	if (s == liteseq::strand::STRAND_FWD)
		return '>';
	else // STRAND_REV
		return '<';
}

std::set<pt::id_t> get_ref_ids_phased(const bd::VG &g, const std::string &sn,
				      pt::u32 phase_idx)
{
	// std::string sn = vcf_file.get_sample_name(sample_idx);
	std::set<pt::id_t> ref_ids = g.get_refs_in_sample(sn);

	std::set<pt::id_t> filtered_ref_ids;
	for (pt::id_t r_id : ref_ids) {
		const povu::refs::Ref &r = g.get_ref_by_id(r_id);

		// TODO: find a better way to handle non PANSN
		if (r.get_format() != PN) // just trust it
			filtered_ref_ids.insert(r_id);
		else if (r.get_format() == PN && r.get_hap_id() == phase_idx)
			filtered_ref_ids.insert(r_id);
	}

	return filtered_ref_ids;
}

std::set<pt::id_t> get_ref_ids(const bd::VG &g,
			       const mto::from_vcf::VCFile &vcf_file,
			       pt::u32 sample_idx, pt::u32 phase_idx)
{
	const std::string &sn = vcf_file.get_sample_name(sample_idx);

	if (g.get_ploidy(sn) == 1 || g.get_ploidy(sn) == pc::INVALID_IDX)
		return g.get_refs_in_sample(sn);
	else { // ploidy is never 0, the else is always >1
		pt::u32 ploidy_id = g.get_ploidy_id(sn, phase_idx);
		return get_ref_ids_phased(g, sn, ploidy_id);
	}
}
} // namespace zien::common
