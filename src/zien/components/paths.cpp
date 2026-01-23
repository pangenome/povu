#include <algorithm>

#include <liteseq/refs.h> // for ref_walk, ref

// #include "povu/common/core.hpp"
#include "povu/common/core.hpp"
#include "povu/graph/bidirected.hpp"	  // for VG
#include "zien/components/components.hpp" // for display_lines

namespace zien::components::paths
{
namespace lq = liteseq;

// TODO: use utils
char lq_strand_to_or_e(lq::strand s)
{
	return (s == lq::strand::STRAND_FWD) ? '>' : '<';
}

void update_paths(const bd::VG &g, display_lines &pd)
{
	std::string curr_l;
	// pt::u32 longest_tag;
	for (pt::u32 h_idx{}; h_idx < g.get_hap_count(); h_idx++) {
		curr_l.clear();
		curr_l += g.get_tag(h_idx);

		if (curr_l.length() > pd.lh)
			pd.lh = curr_l.length();

		pd.meta[h_idx].ref_name_pos = curr_l.length();

		const lq::ref_walk *rw = g.get_ref_vec(h_idx)->walk;

		pt::u32 limit = 1000;
		pt::u32 N = std::min(rw->step_count, limit);

		for (pt::u32 i{}; i < N; i++) {
			pt::id_t v_id = rw->v_ids[i];
			char o = lq_strand_to_or_e(rw->strands[i]);

			curr_l += o;
			curr_l += std::to_string(v_id);
		}
		pd.lines.push_back(curr_l);
	}
}
} // namespace zien::components::paths
