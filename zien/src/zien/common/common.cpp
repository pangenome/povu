#include "zien/common/common.hpp"

#include <set>
#include <string>

#include <liteseq/refs.h>  // for ref_walk, ref
#include <quilt/types.hpp> // for qt

#include "mto/from_vcf.hpp"	     // for VCFile
#include "povu/graph/bidirected.hpp" // for bidirected

#include <povu/refs/refs.hpp> // for Ref

namespace zien::common
{

char to_char(liteseq::strand s)
{
	if (s == liteseq::strand::STRAND_FWD)
		return '>';
	else // STRAND_REV
		return '<';
}

std::set<qt::id_t> get_ref_ids_phased(const bd::VG &g, const std::string &sn,
				      qt::u32 phase_idx)
{
	// std::string sn = vcf_file.get_sample_name(sample_idx);
	std::set<qt::id_t> ref_ids = g.get_refs_in_sample(sn);

	std::set<qt::id_t> filtered_ref_ids;
	for (qt::id_t r_id : ref_ids) {
		const oza::refs::Ref &r = g.get_ref_by_id(r_id);

		// TODO: find a better way to handle non PANSN
		if (r.get_format() != PN) // just trust it
			filtered_ref_ids.insert(r_id);
		else if (r.get_format() == PN && r.get_hap_id() == phase_idx)
			filtered_ref_ids.insert(r_id);
	}

	return filtered_ref_ids;
}

std::set<qt::id_t> get_ref_ids(const bd::VG &g,
			       const mto::from_vcf::VCFile &vcf_file,
			       qt::u32 sample_idx, qt::u32 phase_idx)
{
	const std::string &sn = vcf_file.get_sample_name(sample_idx);

	if (g.get_ploidy(sn) == 1 || g.get_ploidy(sn) == pc::INVALID_IDX)
		return g.get_refs_in_sample(sn);
	else { // ploidy is never 0, the else is always >1
		qt::u32 ploidy_id = g.get_ploidy_id(sn, phase_idx);
		return get_ref_ids_phased(g, sn, ploidy_id);
	}
}
} // namespace zien::common
