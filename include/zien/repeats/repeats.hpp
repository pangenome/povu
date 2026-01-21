#ifndef ZIEN_REPEATS_HPP
#define ZIEN_REPEATS_HPP

#include <liteseq/refs.h> // for ref_walk, ref

// #include "ita/align/align.hpp"	   // for align, aln_level_e
// #include "ita/genomics/allele.hpp" // for ia::at_itn

#include "povu/common/core.hpp" // for pt
// #include "povu/common/utils.hpp"     // for pu::concat_with
#include "povu/graph/bidirected.hpp" // for VG
#include "povu/graph/types.hpp"

namespace zien::tui::repeats
{
namespace lq = liteseq;

std::tuple<std::vector<pt::u32>, std::vector<pt::slice>,
	   std::vector<std::string>>
foo(const bd::VG &g, pt::u32 ref_h_idx, const pt::op_t<ptg::id_or_t> &ef,
    pt::u32 pos);

} // namespace zien::tui::repeats

#endif // ZIEN_STATE_HPP
