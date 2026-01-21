#ifndef ZIEN_COMPONENTS_VCF_HPP
#define ZIEN_COMPONENTS_VCF_HPP

#include <ncurses.h>

#include "mto/from_vcf.hpp"		  // for VCFile
#include "povu/graph/bidirected.hpp"	  // for VG
#include "zien/components/components.hpp" // for Mode

namespace zien::components::vcf
{
/* top left pane */
void comp_vcfs(const bd::VG &g, const mto::from_vcf::VCFile &vcf_file,
	       const std::vector<pt::u32> &invalid_recs, display_lines &pd);

} // namespace zien::components::vcf

#endif // ZIEN_COMPONENTS_VCF_HPP
