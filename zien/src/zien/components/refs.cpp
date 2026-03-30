#include <liteseq/refs.h> // for ref_walk, ref

#include "mto/from_vcf.hpp"		  // for VCFile
#include "povu/graph/bidirected.hpp"	  // for VG
#include "zien/components/common.hpp"	  // for comp_update_refs
#include "zien/components/components.hpp" // for status_bar

namespace zien::components::refs
{

/* top left (ref) pane */
void update_refs(const bd::VG &g, const mto::from_vcf::VCFile &vcf_file,
		 pt::u32 selected_rec_idx, display_lines &pd)
{
	const pt::u32 REF_AT_IDX{0};
	pd.reset();
	zien::components::common::comp_update_refs(
		g, vcf_file, selected_rec_idx, REF_AT_IDX, pd);
};
} // namespace zien::components::refs
