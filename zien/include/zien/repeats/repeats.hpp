#ifndef ZIEN_REPEATS_HPP
#define ZIEN_REPEATS_HPP

#include <vector>

#include <liteseq/refs.h>	    // for ref_walk, ref
#include <oza/graph/bidirected.hpp> // for VG
#include <quilt/graph_types.hpp>    // for v_end_e, side_n_id_t, side_n_idx_t
#include <quilt/types.hpp>	    // for qt

namespace zien::tui::repeats
{
namespace lq = liteseq;

std::tuple<std::vector<qt::u32>, std::vector<qt::slice>,
	   std::vector<std::string>>
foo(const bd::VG &g, qt::u32 ref_h_idx, const qt::op_t<ptg::id_or_t> &ef,
    qt::u32 pos);

} // namespace zien::tui::repeats

#endif // ZIEN_STATE_HPP
