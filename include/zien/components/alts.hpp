#ifndef ZIEN_COMPONENTS_ALTS_HPP
#define ZIEN_COMPONENTS_ALTS_HPP

#include <ncurses.h>

#include "mto/from_vcf.hpp"		  // for VCFile
#include "povu/graph/bidirected.hpp"	  // for VG
#include "zien/components/components.hpp" // for Mode

namespace zien::components::alts
{
/* right (alt) pane */
void update_alts(const bd::VG &g, const mto::from_vcf::VCFile &vcf_file,
		 pt::u32 selected_rec_idx, display_lines &pd);
} // namespace zien::components::alts

#endif // ZIEN_COMPONENTS_ALTS_HPP
