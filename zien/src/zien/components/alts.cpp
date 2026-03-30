#include <ncurses.h>

#include "mto/from_vcf.hpp"		  // for VCFile
#include "povu/graph/bidirected.hpp"	  // for VG
#include "zien/components/common.hpp"	  // for comp_update_refs
#include "zien/components/components.hpp" // for Mode

namespace zien::components::alts
{
/* right (alt) pane */
void update_alts(const bd::VG &g, const mto::from_vcf::VCFile &vcf_file,
		 pt::u32 selected_rec_idx, display_lines &pd)
{
	pd.reset();

	const mto::from_vcf::VCFRecord &rec =
		vcf_file.get_records().at(selected_rec_idx);

	pt::u32 N = rec.get_alt_count() + 1;

	for (pt::u32 i{1}; i < N; i++) { // start at 1 to skip REF
		zien::components::common::comp_update_refs(
			g, vcf_file, selected_rec_idx, i, pd);

		if (i + 1 < N) // add a separator line between alleles
			pd.group_lines.insert(pd.lines.size());
	}
};
} // namespace zien::components::alts
