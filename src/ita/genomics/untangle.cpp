#include "ita/genomics/untangle.hpp"

#include <iterator>
#include <queue>   // for queue
#include <set>	   // for set, operator!=
#include <string>  // for basic_string, string
#include <utility> // for move
#include <vector>  // for vector

#include "ita/align/align.hpp"	   // for align, aln_level_e
#include "ita/genomics/allele.hpp" // for Exp, itn_t

#include "povu/common/core.hpp"	     // for pt, id_t, up_t, operator<
#include "povu/graph/bidirected.hpp" // for VG, bd
#include "povu/graph/types.hpp"	     // for or_e, id_or_t

namespace ita::untangle
{
using namespace ia;

constexpr auto fwd = pgt::or_e::forward;
constexpr auto rev = pgt::or_e::reverse;

pgt::or_e lq_strand_to_or_e(lq::strand s)
{
	return s == lq::strand::STRAND_FWD ? fwd : rev;
}

broad_walk lineup(const bd::VG &g, const std::vector<pt::u32> &sorted_w,
		  pt::u32 h_idx)
{
	broad_walk unrolled;

	const pt::u32 J = sorted_w.size();

	// add after the last position greater than current
	auto add_unrolled = [&](pt::u32 pos, pt::u32 v_id, pgt::or_e o)
	{
		extended_step es{pos, v_id, o};
		auto it = std::upper_bound(
			unrolled.begin(), unrolled.end(), es,
			[](const extended_step &a, const extended_step &b)
			{ return std::get<0>(a) < std::get<0>(b); });
		unrolled.insert(it, es);
	};

	const lq::ref_walk *h_w = g.get_ref_vec(h_idx)->walk; // the hap walk
	for (pt::u32 j{}; j < J; j++) {
		pt::u32 v_id = sorted_w[j];
		const std::vector<pt::u32> &positions =
			g.get_vertex_ref_idxs(g.v_id_to_idx(v_id), h_idx);

		for (pt::u32 i : positions) { // index in the hap walk
			pgt::or_e o = lq_strand_to_or_e(h_w->strands[i]);
			add_unrolled(i, v_id, o);
		}
	}

	return unrolled;
}

race cluster(const broad_walk &unrolled)
{
	/* cluster by the ends of the walk w  */
	race clusters;
	std::queue<broad_step> q;
	broad_walk cluster;
	auto flush_q = [&]()
	{
		if (q.empty())
			return;

		while (!q.empty()) {
			cluster.emplace_back(q.front());
			q.pop();
		}
		clusters.emplace_back(cluster);
		cluster.clear();
	};

	for (auto it = unrolled.begin(); it != unrolled.end(); ++it) {
		q.push(*it);
		if (std::next(it) == unrolled.end()) { // last element
			flush_q();
			break; // TODO[C] return here?
		}

		auto [idx_in_hap_curr, v_id_curr, o] = *it;
		auto [idx_in_hap_nxt, v_id_nxt, _] = *std::next(it);

		if (idx_in_hap_curr + 1 != idx_in_hap_nxt)
			flush_q();
	}

	return clusters;
}

race gen_race(const bd::VG &g, const std::vector<pt::u32> &sorted_w,
	      pt::u32 h_idx)
{
	broad_walk bw = lineup(g, sorted_w, h_idx);
	return cluster(bw);
}

// generate a ia::at_it from a race
// narrow down an itn from a race
ia::at_itn race_to_at_itn(const race &r)
{
	std::vector<pgt::walk_t> allele_traversals;
	pgt::walk_t w;
	for (const ext_lap &lap : r) {
		for (auto &[_, v_id, o] : lap)
			w.emplace_back(pgt::id_or_t{v_id, o});

		allele_traversals.emplace_back(w);
		w.clear();
	}
	ia::at_itn itn(std::move(allele_traversals));

	return itn;
};

void fill_row(const ir::RoV &rov, const race &r, pt::u32 loop_no, pt::u32 h_idx,
	      depth_matrix &dm)
{
	const broad_walk &ref_at = r[loop_no];
	dm.set_loop_no(h_idx, loop_no);
	for (auto &[_, v_id, o] : ref_at) {
		switch (o) {
		case pgt::or_e::forward:
			dm.set_data(h_idx, rov.get_sorted_pos(v_id), 1);
			break;
		case pgt::or_e::reverse:
			dm.set_data(h_idx, rov.get_sorted_pos(v_id), 2);
			break;
		}
	}
}

std::string do_align(const bd::VG &g, const std::vector<pt::u32> &sorted_w,
		     pt::u32 ref_h_idx, pt::u32 alt_h_idx)
{

	race ref_race = gen_race(g, sorted_w, ref_h_idx);
	ia::at_itn ref_itn = race_to_at_itn(ref_race);
	race alt_race = gen_race(g, sorted_w, alt_h_idx);
	ia::at_itn alt_itn = race_to_at_itn(alt_race);

	auto lvl = ita::align::aln_level_e::at;

	std::string et = ita::align::align(ref_itn, alt_itn, lvl);

	return et;
}

bool col_has_mismatch(const std::map<pt::u32, std::string> &alns, pt::u32 j)
{
	for (const auto &[_, aln] : alns)
		if (aln[j] == 'X')
			return true;

	return false;
}

pt::u32 comp_loop_no(pt::u32 loop_col_no, const std::string &et)
{
	pt::u32 loop_no{};
	for (pt::u32 j{}; j < loop_col_no; j++) {
		if (et[j] == 'I')
			continue;

		loop_no++; // increment loop no only on M, I, X
	}

	return loop_no;
}

std::vector<depth_matrix> unroll(const bd::VG &g, const ir::RoV &rov,
				 const std::vector<pt::id_t> &sorted_w,
				 const std::map<pt::u32, std::string> &alns,
				 const race &ref_race, pt::u32 ref_h_idx,
				 const pt::u32 ref_loop_count, pt::u32 I,
				 pt::u32 J)
{
	std::vector<depth_matrix> unrolled_dms;

	for (pt::u32 j{}; j < ref_loop_count; j++) {
		if (col_has_mismatch(alns, j)) {
			// std::cerr << "col " << j << " has X\n";
			depth_matrix dm(I, J);

			fill_row(rov, ref_race, j, ref_h_idx, dm);
			for (const auto &[h_idx, aln] : alns) {
				if (aln[j] == 'D' || aln[j] == 'I')
					continue;

				pt::u32 ln = comp_loop_no(j, aln);
				race r = gen_race(g, sorted_w, h_idx);
				fill_row(rov, r, ln, h_idx, dm);
			}

			dm.set_tangled(true);

			unrolled_dms.emplace_back(std::move(dm));
		}
	}

	return unrolled_dms;
}

std::vector<depth_matrix> untangle(const bd::VG &g,
				   const std::set<pt::u32> &to_call_ref_ids,
				   const depth_matrix &dm, const ir::RoV &rov)
{
	const pt::u32 I = dm.row_count();
	const pt::u32 J = dm.col_count();
	const std::vector<pt::u32> &sorted_w = rov.get_sorted_vertices();

	// ref idx to unrolled dms
	std::map<pt::u32, std::vector<depth_matrix>> ref_to_unrolled_dms;

	for (pt::u32 ref_h_idx : to_call_ref_ids) {

		race ref_race = gen_race(g, sorted_w, ref_h_idx);
		pt::u32 ref_loop_count = ref_race.size();

		std::map<pt::u32, std::string> alns;

		for (pt::u32 h_idx{}; h_idx < I; h_idx++) {
			if (ref_h_idx == h_idx)
				continue;

			alns[h_idx] = do_align(g, sorted_w, ref_h_idx, h_idx);
		}

		ref_to_unrolled_dms[ref_h_idx] =
			unroll(g, rov, sorted_w, alns, ref_race, ref_h_idx,
			       ref_loop_count, I, J);
	}

	for (auto [_, dms] : ref_to_unrolled_dms)
		for (pt::u32 k{}; k < dms.size(); k++)
			dms[k].set_tangled(true);

	// return the first for now
	return ref_to_unrolled_dms.begin()->second;
}

} // namespace ita::untangle
