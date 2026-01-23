#include <ncurses.h>

#include "mto/from_vcf.hpp"		  // for VCFile
#include "povu/graph/bidirected.hpp"	  // for VG
#include "zien/components/components.hpp" // for Mode

namespace zien::components::vcf
{

std::vector<std::string> comp_vcf_lines(const mto::from_vcf::VCFile &vcf_file)
{
	std::vector<std::string> lines;
	std::stringstream ss; // string buffer

	// header
	const std::vector<std::string> COL_NAMES = {
		"VAR TYPE", "TANGLED", "POS",	 "ID",
		"REF",	    "ALT",     "REF AT", "ALT AT"};

	lines.emplace_back(pu::concat_with(COL_NAMES, '\t'));

	const std::vector<mto::from_vcf::VCFRecord> &all_recs =
		vcf_file.get_records();
	for (const mto::from_vcf::VCFRecord &rec : all_recs) {
		rec.viewer(ss);
		lines.push_back(ss.str());
		ss.str(std::string()); // Clear the stringstream
	}

	return lines;
}

void comp_vcfs(const bd::VG &g, const mto::from_vcf::VCFile &vcf_file,
	       const std::vector<pt::u32> &invalid_recs, display_lines &pd)
{
	/* top left pane */

	std::vector<std::string> vcf_lines = comp_vcf_lines(vcf_file);
	pd.lines.swap(vcf_lines);
	for (auto r_idx : invalid_recs)
		pd.special_lines.insert(r_idx + 1); // +1 for header line
}

} // namespace zien::components::vcf
