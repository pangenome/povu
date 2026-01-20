#ifndef ZIEN_COMPONENTS_GT_HPP
#define ZIEN_COMPONENTS_GT_HPP

#include <ncurses.h>

#include "mto/from_vcf.hpp"		  // for VCFile
#include "povu/graph/bidirected.hpp"	  // for VG
#include "zien/components/components.hpp" // for Mode

namespace zien::components::genotypes
{
/* top-right (haplotype) pane */
void update_haps(const mto::from_vcf::VCFile &vcf_file,
		 pt::u32 selected_rec_idx, display_lines &pd);

} // namespace zien::components::genotypes

#endif // ZIEN_COMPONENTS_HPP
