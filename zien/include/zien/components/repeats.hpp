#ifndef ZIEN_COMPONENTS_REPEATS_HPP
#define ZIEN_COMPONENTS_REPEATS_HPP

#include <ncurses.h>

#include <mto/from_vcf.hpp>		  // for VCFile
#include <oza/graph/bidirected.hpp>	  // for VG
#include <quilt/types.hpp>		  // for qt
					  //
#include "zien/components/components.hpp" // for Mode

namespace zien::components::repeats
{
void update_repeats(const bd::VG &g, const mto::from_vcf::VCFile &vcf_file,
		    qt::u32 rec_idx, display_lines &pd);
} // namespace zien::components::repeats

#endif // ZIEN_COMPONENTS_REPEATS_HPP
