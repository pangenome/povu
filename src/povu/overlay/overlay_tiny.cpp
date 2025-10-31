#include "povu/overlay/overlay.hpp"

#include <utility> // for pair

#include <liteseq/refs.h> // for ref_walk, ref

#include "povu/common/constants.hpp"
#include "povu/common/core.hpp"
#include "povu/genomics/allele.hpp"
#include "povu/graph/types.hpp"
#include "povu/variation/rov.hpp"

namespace povu::overlay
{
namespace lq = liteseq;

constexpr pvr::var_type_e sub = pvr::var_type_e::sub;
constexpr pgt::or_e fo = pgt::or_e::forward;
constexpr pgt::or_e ro = pgt::or_e::reverse;

bool extend_right(const lq::ref_walk *ref_w, const pgt::walk_t &graph_w,
		  pt::u32 graph_w_start_idx, pt::u32 ref_w_start_idx,
		  pt::u32 len)
{

	pt::u32 i = ref_w_start_idx;
	pt::u32 j = graph_w_start_idx;

	pt::u32 g_end = graph_w_start_idx + len;
	for (; j < g_end; i++, j++) {
		pt::idx_t ref_v_id = ref_w->v_ids[i];
		pgt::or_e ref_o = ref_w->strands[i] == lq::strand::STRAND_FWD
					  ? pgt::or_e::forward
					  : pgt::or_e::reverse;

		auto [g_v_id, g_o] = graph_w[j];

		if (ref_v_id != g_v_id || ref_o != g_o)
			return false;
	}

	return true;
}

bool extend_left(const lq::ref_walk *ref_w, const pgt::walk_t &graph_w,
		 pt::u32 graph_w_start_idx, pt::u32 ref_w_start_idx,
		 pt::u32 len)
{
	pt::u32 i = ref_w_start_idx; // go backwards
	pt::u32 j = graph_w_start_idx;

	pt::u32 g_end = graph_w_start_idx + len;
	for (; j < g_end && i != pc::INVALID_IDX; i--, j++) {

		pt::idx_t ref_v_id = ref_w->v_ids[i];
		pgt::or_e ref_o = ref_w->strands[i] == lq::strand::STRAND_FWD
					  ? pgt::or_e::forward
					  : pgt::or_e::reverse;

		auto g_v_id = graph_w[j].v_id;
		auto g_o = pgt::flip(graph_w[j].orientation);

		if (ref_v_id != g_v_id || ref_o != g_o)
			return false;
	}

	if (j == g_end)
		return true;

	return false;
}

void update_exp(pt::u32 r_idx, pt::u32 w_idx, pga::allele_slice_t &&at,
		pga::Exp &e)
{
	std::map<pt::id_t, pga::itn_t> &ref_map = e.get_ref_itns_mut();
	pga::itn_t &itn = ref_map[r_idx];

	std::map<pt::idx_t, std::set<pt::idx_t>> &walk_to_refs =
		e.get_walk_to_ref_idxs_mut();

	itn.append_at_sorted(std::move(at));
	walk_to_refs[w_idx].insert(r_idx);
}

std::pair<std::vector<pga::Exp>, std::vector<sub_inv>>
overlay_tiny(const bd::VG &g, const pvr::RoV &rov,
	     const std::set<pt::id_t> &to_call_ref_ids)
{
	pga::Exp e(&rov);

	bool dbg = (rov.as_str() == ">188>190") ? true : false;

	// print walks
	// if (dbg) {
	//	std::cerr << "RoV Walks:\n";
	//	for (auto &w : rov.get_walks()) {
	//		std::cerr << ptg::to_string(w) << "\n";
	//	}
	//	std::exit(EXIT_FAILURE);
	// }

	const std::vector<pgt::walk_t> &walks = rov.get_walks();
	const pvr::var_type_e DEFAULT_VT = sub;
	const pt::u32 WALK_COUNT = walks.size();
	const pt::u32 WALK_START{};
	// number of times a ref traverses the RoV
	std::map<pt::u32, pt::u32> ref_loop_count;

	for (pt::u32 w_idx{}; w_idx < WALK_COUNT; ++w_idx) {
		const pgt::walk_t &w = walks[w_idx];
		auto [v_id, o] = w.front();

		const pt::u32 w_len = w.size();
		const std::vector<std::vector<pt::idx_t>> &vtx_ref_idxs =
			g.get_vertex_refs(v_id);

		for (pt::u32 r_idx{}; r_idx < g.ref_count(); r_idx++) {
			const std::vector<pt::u32> &r_starts =
				vtx_ref_idxs[r_idx];

			if (r_starts.empty())
				continue;

			const lq::ref_walk *ref_w = g.get_ref_vec(r_idx)->walk;

			for (pt::u32 ref_start : r_starts) {

				bool is_right = extend_right(
					ref_w, w, WALK_START, ref_start, w_len);

				bool is_left = extend_left(ref_w, w, WALK_START,
							   ref_start, w_len);

				if (!is_right && !is_left)
					continue;

				pgt::or_e at_or = is_right ? fo : ro;

				pga::allele_slice_t at{
					&w,    w_idx, WALK_START,
					ref_w, r_idx, ref_start,
					w_len, at_or, DEFAULT_VT};

				update_exp(r_idx, w_idx, std::move(at), e);

				ref_loop_count[r_idx]++;
			}
		}
	}

	for (auto &[_, count] : ref_loop_count)
		if (count > 1)
			e.set_tangled(true);

	return {{e}, {}};
}

} // namespace povu::overlay
