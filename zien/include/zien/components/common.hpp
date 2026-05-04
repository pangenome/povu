#ifndef ZIEN_COMPONENTS_REFS_HPP
#define ZIEN_COMPONENTS_REFS_HPP

#include <ncurses.h>

#include <mto/from_vcf.hpp>	    // for VCFile
#include <oza/graph/bidirected.hpp> // for VG
#include <quilt/types.hpp>	    // for qt

#include "zien/components/components.hpp" // for Mode

namespace zien::components::common
{
void comp_update_refs(const bd::VG &g, const mto::from_vcf::VCFile &vcf_file,
		      qt::u32 selected_rec_idx, qt::u32 at_idx,
		      display_lines &pd);

} // namespace zien::components::common

#endif // ZIEN_COMPONENTS_HPP
