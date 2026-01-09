#include "povu/overlay/overlay.hpp" // for pos, pin_cushion

#include <algorithm>
#include <cstdlib> // for exit, EXIT_FAILURE
#include <map>	   // for map
#include <optional>
#include <ostream>
#include <queue> // for queue
#include <set>
#include <tuple>
#include <utility>
#include <vector> // for vector

#include <liteseq/refs.h> // for ref_walk, ref

#include "povu/common/constants.hpp"
#include "povu/common/core.hpp"
#include "povu/common/log.hpp"
#include "povu/common/utils.hpp"
#include "povu/genomics/allele.hpp"
#include "povu/genomics/untangle.hpp"
#include "povu/graph/types.hpp"
#include "povu/overlay/sne.hpp"
#include "povu/variation/rov.hpp"

namespace povu::overlay
{

using namespace povu::genomics::allele;

bool is_row_empty(const std::vector<pt::u32> &row)
{
	for (pt::u32 e : row)
		if (e > 0)
			return false;

	return true;
}

using bounds_to_alts = std::map<pt::op_t<pt::u32>, std::set<pt::u32>>;

std::vector<pt::op_t<pt::u32>>
find_pairwise_rov(const std::vector<pt::u32> &ref_row,
		  const std::vector<pt::u32> &alt_row, const pt::u32 J)
{
	if (J == 1)
		return {};

	std::vector<pt::op_t<pt::u32>> bounds;
	std::queue<pt::u32> q; // a buffer to store similar cols
	pt::u32 u{pc::INVALID_IDX}, v{pc::INVALID_IDX};
	for (pt::u32 j{}; j < J; j++) {
		if (ref_row[j] == alt_row[j] && ref_row[j] > 0)
			q.push(j);

		if (q.size() > 1) {
			u = q.front();
			q.pop();
			v = q.front();

			// if (u,v) are not adjacent & contain a difference
			for (pt::u32 i{u + 1}; i < v; i++)
				if (ref_row[i] != alt_row[i]) {
					bounds.emplace_back(u, v);
					break;
				}
		}
	}

	return bounds;
}

bool match_in_cxt(const std::vector<pt::u32> &ref_row,
		  const std::vector<pt::u32> &alt_row, const pt::u32 j_u,
		  const pt::u32 j_v)
{
	for (pt::u32 j{j_u}; j <= j_v; j++)
		if (ref_row[j] != alt_row[j])
			return false;

	return true;
}

bool is_inverted(const std::vector<pt::u32> &ref_row,
		 const std::vector<pt::u32> &alt_row, const pt::u32 J)
{
	if (J == 1)
		return false;

	// returns true if col j is inverted
	auto is_col_inverted = [&](pt::u32 j) -> bool
	{
		return (ref_row[j] == 1 && alt_row[j] == 2) ||
		       (ref_row[j] == 2 && alt_row[j] == 1);
	};

	for (pt::u32 j{}; j < J; j++)
		if (!is_col_inverted(j))
			return false;

	return true;
}

std::vector<pt::op_t<pt::u32>> comp_bounds(const std::vector<pt::u32> &ref_row,
					   const std::vector<pt::u32> &alt_row,
					   const pt::u32 J, bool dbg)
{
	if (J == 1)
		return {};

	auto comp_bounds = [&](pt::u32 j) -> bool
	{
		return ref_row[j] > 0 && (ref_row[j] == alt_row[j] ||
					  ref_row[j] == 1 && alt_row[j] == 2 ||
					  ref_row[j] == 2 && alt_row[j] == 1);
	};

	std::vector<pt::op_t<pt::u32>> bounds;
	std::queue<pt::u32> q; // a buffer to store similar cols
	pt::u32 u{pc::INVALID_IDX}, v{pc::INVALID_IDX};
	for (pt::u32 j{}; j < J; j++) {
		// if (ref_row[j] == alt_row[j] && ref_row[j] > 0)
		//	q.push(j);

		if (comp_bounds(j))
			q.push(j);

		if (q.size() > 1) {
			u = q.front();
			q.pop();
			v = q.front();

			// if (u,v) are not adjacent & contain a difference
			for (pt::u32 i{u + 1}; i < v; i++)
				if (ref_row[i] != alt_row[i]) {
					bounds.emplace_back(u, v);
					break;
				}
		}
	}

	return bounds;
}

bool slice_match(const pga::hap_slice &a, const pga::hap_slice &b)
{
	if (a.len != b.len)
		return false;

	for (pt::u32 i{}; i < a.len; i++) {
		auto step_a = a.get_step(a.ref_start_idx + i);
		auto step_b = b.get_step(b.ref_start_idx + i);

		if (step_a != step_b)
			return false;
	}

	return true;
}

bool is_inv_local(const pga::hap_slice &a, const pga::hap_slice &b)
{
	if (a.len != b.len)
		return false;

	for (pt::u32 i{}, j{b.len - 1}; i < a.len; i++, j--) {
		auto step_a = a.get_step(a.ref_start_idx + i);
		auto step_b = b.get_step(b.ref_start_idx + j);

		// invert strand
		auto [v_id_b, o_b] = step_b;
		pgt::or_e o_b_inv = o_b == pgt::or_e::forward
					    ? pgt::or_e::reverse
					    : pgt::or_e::forward;

		if (step_a != pgt::step_t{v_id_b, o_b_inv})
			return false;
	}

	return true;
}

pga::hap_slice hap_sl_from_lap(const bd::VG &g, const pga::rov_boundaries &cxt,
			       const ext_lap &lap, pt::u32 h_idx, bool dbg)
{
	auto [u, v] = cxt.get_bounds();
	auto [l_cxt_v_id, _] = u;
	auto [r_cxt_v_id, __] = v;

	std::vector<pt::u32> starts;
	std::vector<pt::u32> ends;

	auto update_q = [](pt::u32 idx_in_hap, std::vector<pt::u32> &q)
	{
		// insert in sorted order by min
		auto it = std::lower_bound(q.begin(), q.end(), idx_in_hap);
		q.insert(it, idx_in_hap);
	};

	// std::optional<pt::u32> start_opt, end_opt;
	for (const auto &[idx_in_hap, v_id, ___] : lap) {
		if (l_cxt_v_id == v_id)
			update_q(idx_in_hap, starts);

		// start_opt = idx_in_hap;

		if (r_cxt_v_id == v_id)
			update_q(idx_in_hap, ends);
		// end_opt = idx_in_hap;

		// if (start_opt && end_opt) // done
		//	break;
	}

	if (starts.empty() || ends.empty()) {
		ERR("Could not find pga::context nodes in haplotype {} {}",
		    h_idx, cxt.to_string());
		std::exit(EXIT_FAILURE);
	}

	if (starts.front() > ends.front()) // swap the starts and ends
		std::swap(starts, ends);

	pt::u32 start = starts.front();
	pt::u32 end = ends.back();
	// pt::u32 start = std::min(*start_opt, *end_opt);
	// pt::u32 end = std::max(*start_opt, *end_opt);
	pt::u32 len = end - start + 1;

	if (dbg) {
		std::cerr << "starts " << pu::concat_with(starts, ',') << "\n";
		std::cerr << "ends " << pu::concat_with(ends, ',') << "\n";
		INFO("hap_sl_from_lap for haplotype {} cxt {}: start "
		     "{} end {} "
		     "len {}",
		     h_idx, cxt.to_string(), start, end, len);
	}

	return pga::hap_slice{g.get_ref_vec(h_idx)->walk, h_idx, start, len};
}

std::optional<pga::hap_slice>
hap_sl_from_lap_alt(const bd::VG &g, const pga::rov_boundaries &cxt,
		    const ext_lap &lap, pt::u32 h_idx, bool dbg)
{
	auto [u, v] = cxt.get_bounds();
	auto [l_cxt_v_id, _] = u;
	auto [r_cxt_v_id, __] = v;

	std::vector<pt::u32> starts;
	std::vector<pt::u32> ends;

	auto update_q = [](pt::u32 idx_in_hap, std::vector<pt::u32> &q)
	{
		// insert in sorted order by min
		auto it = std::lower_bound(q.begin(), q.end(), idx_in_hap);
		q.insert(it, idx_in_hap);
	};

	// std::optional<pt::u32> start_opt, end_opt;
	for (const auto &[idx_in_hap, v_id, ___] : lap) {
		if (l_cxt_v_id == v_id)
			update_q(idx_in_hap, starts);

		// start_opt = idx_in_hap;

		if (r_cxt_v_id == v_id)
			update_q(idx_in_hap, ends);
		// end_opt = idx_in_hap;

		// if (start_opt && end_opt) // done
		//	break;
	}

	if (starts.empty() || ends.empty()) {
		return std::nullopt;
	}

	if (starts.front() > ends.front()) // swap the starts and ends
		std::swap(starts, ends);

	pt::u32 start = starts.front();
	pt::u32 end = ends.back();
	// pt::u32 start = std::min(*start_opt, *end_opt);
	// pt::u32 end = std::max(*start_opt, *end_opt);
	pt::u32 len = end - start + 1;

	if (dbg) {
		std::cerr << "starts " << pu::concat_with(starts, ',') << "\n";
		std::cerr << "ends " << pu::concat_with(ends, ',') << "\n";
		INFO("hap_sl_from_lap for haplotype {} cxt {}: start "
		     "{} end {} "
		     "len {}",
		     h_idx, cxt.to_string(), start, end, len);
	}

	return pga::hap_slice{g.get_ref_vec(h_idx)->walk, h_idx, start, len};
}

pga::rov_boundaries gen_cxt(pt::u32 u, pt::u32 v, const ext_lap &lap)

{
	std::vector<extended_step> u_buff;
	std::vector<extended_step> v_buff;

	for (pt::u32 i{}; i < lap.size(); i++) {
		auto [_, v_id, __] = lap[i];
		if (v_id == u)
			u_buff.push_back(lap[i]);

		if (v_id == v)
			v_buff.push_back(lap[i]);
	}

	auto [_, l_v_id, l_o] = u_buff.front();
	auto [__, r_v_id, r_o] = v_buff.back();

	return pga::rov_boundaries{pgt::id_or_t{l_v_id, l_o},
				   pgt::id_or_t{r_v_id, r_o}};
}

pos::pin gen_pin(pt::u32 h_idx, const extended_step &s)
{
	auto [idx_in_hap, _, __] = s;

	return {h_idx, idx_in_hap};
}

void print_race(std::ostream &os, const race &r)
{
	for (const auto &lap : r) {
		for (auto [_, v_id, o] : lap) {
			pgt::id_or_t v{v_id, o};
			std::cerr << v;
		}
		os << "\t";
	}
	os << "\n";
}

pga::trek comp_exps(const bd::VG &g, const pvr::RoV *rov,
		    const std::set<pt::u32> &to_call_ref_ids,
		    const depth_matrix &dm, bool tangled,
		    pos::pin_cushion &pcushion, bool dbg)
{
	const pt::u32 I = dm.row_count();
	const pt::u32 J = rov->size();

	auto tk = pga::trek::create_new(rov, I, tangled);

	// memoise race computation
	std::map<pt::u32, race> h_idx_to_race_map;
	auto h_idx_to_race = [&](pt::u32 h_idx) -> const race &
	{
		if (pv_cmp::contains(h_idx_to_race_map, h_idx))
			return h_idx_to_race_map.at(h_idx);

		const std::vector<pt::id_t> &sorted_w =
			rov->get_sorted_vertices();

		h_idx_to_race_map[h_idx] =
			povu::genomics::untangle::gen_race(g, sorted_w, h_idx);

		// print race
		if (dbg) {
			std::cerr << __func__ << " race for haps " << h_idx
				  << ":\n";
			for (const auto &lap : h_idx_to_race_map.at(h_idx)) {
				for (const auto &[_, v_id, o] : lap) {
					pgt::id_or_t v{v_id, o};
					std::cerr << v << ",";
				}
				std::cerr << "\t";
			}
			std::cerr << "\n";
		}

		return h_idx_to_race_map.at(h_idx);
	};

	// auto init_context = [&](pt::u32 ref_h_idx,
	//			const pga::rov_boundaries cxt,
	//			const ext_lap &ref_lap)
	// {
	//	std::map<pga::rov_boundaries, pga::minimal_rov> &d =
	//		tk.get_ref_recs_mut(ref_h_idx);

	//	d.insert({cxt,
	//		  {cxt, hap_sl_from_lap(g, cxt, ref_lap,
	// ref_h_idx)}});
	// };

	// auto add_alt = [&](pt::u32 ref_h_idx, pt::u32 alt_h_idx,
	//		   const pga::rov_boundaries cxt,
	//		   const ext_lap &alt_lap)
	// {
	//	pga::minimal_rov &min_rov =
	//		tk.get_ref_recs_mut(ref_h_idx).at(cxt);

	//	min_rov.add_alt(hap_sl_from_lap(g, cxt, alt_lap,
	// alt_h_idx));
	// };

	for (auto ref_h_idx : to_call_ref_ids) {
		std::vector<pt::u32> ref_row = dm.get_row_data(ref_h_idx);

		for (pt::u32 h_idx{}; h_idx < dm.row_count(); h_idx++) {
			if (ref_h_idx == h_idx)
				continue;

			std::vector<pt::u32> alt_row = dm.get_row_data(h_idx);
			if (is_row_empty(alt_row)) {
				tk.add_no_cov(ref_h_idx, h_idx);
				continue;
			}

			std::vector<pt::op_t<pt::u32>> cxts =
				comp_bounds(ref_row, alt_row, J, dbg);

			bool is_inv = is_inverted(ref_row, alt_row, J);

			if (dbg) {
				INFO("ref {} alt {}", ref_h_idx, h_idx);
				INFO("{}", is_inv);

				// cxts
				for (auto [u, v] : cxts)
					INFO("  cxt: ({} , {})", u, v);
			}

			if (cxts.empty() && !is_inv) {
				tk.add_match_ref(ref_h_idx, h_idx);
				continue;
			}

			const race &ref_race = h_idx_to_race(ref_h_idx);
			const race &alt_race = h_idx_to_race(h_idx);

			// does nothing if the ref idx is already
			// initialised tk.init_ref_idx(ref_h_idx);

			const ext_lap &ref_lap =
				ref_race[dm.get_loop_no(ref_h_idx)];

			const ext_lap &alt_lap =
				alt_race[dm.get_loop_no(h_idx)];

			if (is_inv) {
				pcushion.add_pin_pair(
					{gen_pin(ref_h_idx, ref_lap.front()),
					 gen_pin(h_idx, alt_lap.back())});
				continue;
			}

			// if (dbg)
			//	volatile int pause = 1;

			// std::cerr << __func__
			//	  << " found contexts for ref haplotype
			//"
			//	  << ref_h_idx << " and alt haplotype "
			//<< h_idx
			//	  << ":\n";

			for (auto b : cxts) {
				pt::id_t u = rov->get_sorted_vertex(b.first);
				pt::id_t v = rov->get_sorted_vertex(b.second);

				pga::rov_boundaries c = gen_cxt(u, v, ref_lap);

				pga::cxt_to_min_rov_map &m =
					tk.get_min_rov(ref_h_idx);

				pga::hap_slice ref_sl = hap_sl_from_lap(
					g, c, ref_lap, ref_h_idx, dbg);

				pga::hap_slice alt_sl = hap_sl_from_lap(
					g, c, alt_lap, h_idx, dbg);

				// TODO: check locally
				if (slice_match(ref_sl, alt_sl) ||
				    is_inv_local(ref_sl, alt_sl)) {
					continue;
				}

				if (!pv_cmp::contains(m, c)) { // init
					// pga::hap_slice ref_sl =
					// hap_sl_from_lap(	g, c, ref_lap,
					// ref_h_idx, dbg);

					m.insert({c,
						  pga::minimal_rov(
							  c,
							  std::move(ref_sl))});
					// m[c] = pga::minimal_rov(
					//	c,
					//	hap_sl_from_lap(g, c,
					// ref_lap,
					// ref_h_idx));
				}

				pga::minimal_rov &min_rov = m.at(c);
				// pga::hap_slice alt_sl = hap_sl_from_lap(
				//	g, c, alt_lap, h_idx, dbg);
				min_rov.add_alt(std::move(alt_sl));
				// if (!tk.has_context(ref_h_idx, c))
				//	init_context(ref_h_idx, c,
				// ref_lap);

				// add_alt(ref_h_idx, h_idx, c,
				// alt_lap);
			}
		}

		if (!tk.has_data())
			return tk;

		for (auto &[cxt, min_rov] : tk.get_ref_recs_mut(ref_h_idx)) {
			auto [l, r] = cxt.get_bounds();
			auto [v_id_l, _] = l;
			auto [v_id_r, __] = r;
			pt::u32 j_u_ = rov->get_sorted_pos(v_id_l);
			pt::u32 j_v_ = rov->get_sorted_pos(v_id_r);

			pt::u32 j_u = std::min(j_u_, j_v_);
			pt::u32 j_v = std::max(j_u_, j_v_);

			const race &ref_race = h_idx_to_race(ref_h_idx);
			const ext_lap &ref_lap =
				ref_race[dm.get_loop_no(ref_h_idx)];
			pga::hap_slice ref_sl = hap_sl_from_lap(g, cxt, ref_lap,
								ref_h_idx, dbg);

			for (pt::u32 h_idx{}; h_idx < I; h_idx++) {
				std::vector<pt::u32> alt_row =
					dm.get_row_data(h_idx);

				// check context match, cheap filter
				if (!match_in_cxt(ref_row, alt_row, j_u, j_v))
					continue;

				const race &alt_race = h_idx_to_race(h_idx);

				if (!(dm.get_loop_no(h_idx) < alt_race.size()))
					continue;

				// std::cerr << alt_race.size() << " "
				//	  << dm.get_loop_no(h_idx) << "\n";

				const ext_lap &alt_lap =
					alt_race[dm.get_loop_no(h_idx)];
				std::optional<pga::hap_slice> opt_alt_sl =
					hap_sl_from_lap_alt(g, cxt, alt_lap,
							    h_idx, dbg);

				if (!opt_alt_sl)
					continue;

				pga::hap_slice alt_sl = opt_alt_sl.value();

				if (slice_match(ref_sl, alt_sl))
					min_rov.add_haps_match_ref(h_idx);
			}
		}
	}

	if (dbg) {
		tk.dbg_print(std::cerr);
	}

	return tk;
}

/**
 * [out] rov_exps: vector of expeditions, one per pairwise variant set
 */
std::vector<pga::trek> overlay_generic(const bd::VG &g, pvr::RoV &rov,
				       const std::set<pt::u32> &to_call_ref_ids,
				       pos::pin_cushion &pcushion)
{
	bool dbg = rov.as_str() == ">1>4" ? true : false;
	if (dbg)
		INFO("{}", rov.as_str());

	std::vector<pga::trek> treks;

	pt::u32 I = g.get_hap_count();	    // rows
	pt::u32 J = rov.get_vertex_count(); // cols

	const std::vector<pt::u32> &sorted_vertices = rov.get_sorted_vertices();

	depth_matrix dm(I, J, sorted_vertices);
	dm.fill(g, rov);

	if (dbg)
		dm.print(std::cerr);

	if (dm.tangled()) {
		std::vector<depth_matrix> unrolled_dms =
			povu::genomics::untangle::untangle(g, to_call_ref_ids,
							   dm, rov);

		for (pt::u32 k{}; k < unrolled_dms.size(); k++)
			unrolled_dms[k].set_tangled(true);

		for (pt::u32 i{}; i < unrolled_dms.size(); i++) {
			const depth_matrix &dm_ = unrolled_dms.at(i);
			pga::trek tk = comp_exps(g, &rov, to_call_ref_ids, dm_,
						 true, pcushion, dbg);
			treks.emplace_back(std::move(tk));
		}
	}
	else {
		pga::trek tk = comp_exps(g, &rov, to_call_ref_ids, dm, false,
					 pcushion, dbg);

		treks.emplace_back(std::move(tk));
	}

	return treks;
}
} // namespace povu::overlay
