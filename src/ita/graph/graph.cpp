#include "ita/graph/graph.hpp"

#include <cstdlib> // for exit
#include <set>	   // for set

#include "ita/graph/bfs_tree.hpp"

#include "povu/common/compat.hpp" // for pv_cmp, contains, format
#include "povu/common/core.hpp"
#include "povu/common/log.hpp" // for WARN, ERR
#include "povu/common/utils.hpp"
#include "povu/graph/types.hpp"

namespace povu::genomics::graph
{

/**
 *@brief get the edges of a vertex end
 *
 *@param v the vertex
 *@param ve the vertex end
 *@return a set of edge indices
 */
inline const std::set<pt::idx_t> edges_at_end(const bd::Vertex &v,
					      pgt::v_end_e ve) noexcept
{
	return ve == pgt::v_end_e::l ? v.get_edges_l() : v.get_edges_r();
}

/**
 *@brief get vertex end from the traversal orientation and direction
 *
 *@param o the orientation
 *@param e the direction
 *@return the end of the vertex
 */
inline pgt::v_end_e get_v_end(pgt::or_e o, dir_e e) noexcept
{
	switch (e) {
	case dir_e::in:
		return (o == pgt::or_e::forward) ? pgt::v_end_e::l
						 : pgt::v_end_e::r;
	case dir_e::out:
		return (o == pgt::or_e::forward) ? pgt::v_end_e::r
						 : pgt::v_end_e::l;
	}

	PL_ERR("Invalid v_end : {}", static_cast<int>(e));
	std::exit(EXIT_FAILURE);
};

/**
 *@brief get the orientation of a vertex end based on the side and direction of
 *traversal
 *
 *@param side the side of the vertex end
 *@param d the direction of traversal
 *@return the orientation of the vertex end
 */
inline pgt::or_e get_or(pgt::v_end_e side, dir_e d) noexcept
{
	switch (d) {
	case IN:
		return (side == pgt::v_end_e::l ? pgt::or_e::forward
						: pgt::or_e::reverse);
	case OUT:
		return (side == pgt::v_end_e::l ? pgt::or_e::reverse
						: pgt::or_e::forward);
	}

	PL_ERR("Invalid orientation: {}", static_cast<int>(d));
	std::exit(EXIT_FAILURE);
};

ita::bfs::BfsTree comp_bfs_tree(const bd::VG &g, pvst::route_e route,
				idx_or_t src, idx_or_t snk)
{
	// default is source to sink
	dir_e ve_dir = OUT;
	dir_e nbr_dir = IN;
	idx_or_t start = src;
	idx_or_t end = snk;

	if (__builtin_expect(route == pvst::route_e::e2s, 0)) {
		ve_dir = IN;
		nbr_dir = OUT;
		start = snk;
		end = src;
	}

	std::queue<idx_or_t> q;
	q.push(start);

	std::set<pt::u32> seen;
	seen.insert(start.v_id);

	ita::bfs::BfsTree t;
	pt::u32 t_v_idx =
		t.add_vertex({g.v_idx_to_id(start.v_id), start.orientation});

	std::map<idx_or_t, pt::u32> parent_map;
	parent_map[start] = t_v_idx;

	t.set_start(t_v_idx);

	while (!q.empty()) {
		// get the incoming vertices based on orientation
		idx_or_t curr = q.front();
		q.pop();

		if (q.size() > MAX_UNBLOCK_CTR) // pressure valve
			continue;

		if (curr == end) {
			t.set_end(parent_map[curr]);
			// continue;
		}

		auto [v_idx, o] = curr;

		pgt::v_end_e ve = get_v_end(o, ve_dir);
		const bd::Vertex &v = g.get_vertex_by_idx(v_idx);
		const std::set<pt::idx_t> &nbr_edges = edges_at_end(v, ve);

		for (pt::u32 e_idx : nbr_edges) {
			const bd::Edge &e = g.get_edge(e_idx);
			auto [side, alt_idx] = e.get_other_vtx(v_idx, ve);
			idx_or_t nbr{alt_idx, get_or(side, nbr_dir)};

			if (pv_cmp::contains(seen, nbr.v_id)) { // cross edge
				pt::u32 from = parent_map[curr];
				pt::u32 to = parent_map[nbr];

				t.add_cross_edge(from, to);
				continue;
			}

			// if (seen.find(nbr.v_id) != seen.end()) {
			//	// cross edge
			//	pt::u32 from = parent_map[curr];
			//	pt::u32 to = parent_map[nbr];

			//	t.add_cross_edge(from, to);
			//	continue;
			// }

			if (curr == end)
				continue;

			pt::u32 t_c_v_idx = t.add_vertex(
				{g.v_idx_to_id(nbr.v_id), nbr.orientation});

			parent_map[nbr] = t_c_v_idx;

			t.add_tree_edge(parent_map[curr], t_c_v_idx);

			seen.insert(nbr.v_id);
			q.push(nbr);
		}
	}

	return t;
}

void find_laps(const bd::VG &g, pt::u32 h_idx, pt::id_t u, pt::id_t v,
	       std::vector<pt::slice> &laps)
{
	std::vector<pt::u32> u_hap_idxs =
		g.get_vertex_ref_idxs(g.v_id_to_idx(u), h_idx);

	std::vector<pt::u32> v_hap_idxs =
		g.get_vertex_ref_idxs(g.v_id_to_idx(v), h_idx);

	if (u_hap_idxs.empty() || v_hap_idxs.empty())
		return;

	// sort
	std::sort(u_hap_idxs.begin(), u_hap_idxs.end());
	std::sort(v_hap_idxs.begin(), v_hap_idxs.end());

	pt::u32 N = std::min((u_hap_idxs.size()), v_hap_idxs.size());

	for (pt::u32 i{}; i < N; i++) {
		pt::u32 start = u_hap_idxs[i];
		pt::u32 end = v_hap_idxs[i];

		if (start > end)
			std::swap(start, end);

		laps.emplace_back(start, end - start + 1);
	}
}

std::list<pt::u32> gen_sort(const bd::VG &g, ir::RoV &rov,
			    const pt::u32 HAP_COUNT, bool dbg)
{
	INFO("Generating sort for RoV: {}", rov.as_str());
	dbg = rov.as_str() == ">96686>96691" ? true : false;
	dbg = false;

	auto [l, r, route] = *rov.get_pvst_vtx()->get_route_params();
	auto [start_id, _] = l;
	auto [stop_id, __] = r;

	std::list<pt::u32> sw;
	std::map<pt::u32, std::list<pt::u32>::iterator> w_to_it;

	std::vector<pt::slice> laps;
	laps.reserve(1024);

	// print each lap for each haplotype
	if (dbg) {
		for (pt::u32 h_idx{}; h_idx < HAP_COUNT; h_idx++) {
			laps.clear();
			find_laps(g, h_idx, start_id, stop_id, laps);

			if (dbg && !laps.empty())
				std::cerr << "\n" << h_idx << "\t";

			const liteseq::ref_walk *rw =
				g.get_ref_vec(h_idx)->walk;

			for (const auto &lap : laps) {
				auto [start, len] = lap.data();

				for (pt::u32 j{start}; j < (start + len); j++) {
					pt::u32 v_id = rw->v_ids[j];

					if (dbg)
						std::cerr << v_id << ",";
				}
			}
		}

		std::cerr << "\n---------------------\n";
	}

	auto comp_steps_cxt =
		[&](const pt::slice &lap, const liteseq::ref_walk *rw,
		    pt::op_t<std::map<pt::u32, std::list<pt::u32>::iterator>>
			    &sw_cxt) -> void
	{
		auto &[left_cxt, right_cxt] = sw_cxt;

		std::optional<std::list<pt::u32>::iterator> x;
		std::optional<std::list<pt::u32>::iterator> y;

		auto [start, len] = lap.data();
		for (pt::u32 j{start}; j < (start + len); j++) {
			pt::u32 v_id = rw->v_ids[j];

			if (pv_cmp::contains(left_cxt, v_id))
				continue;

			if (pv_cmp::contains(w_to_it, v_id))
				left_cxt[v_id] = w_to_it.at(v_id);
			else if (j > start) {
				if (pv_cmp::contains(w_to_it,
						     rw->v_ids[j - 1])) {
					x = w_to_it.at(rw->v_ids[j - 1]);
				}

				left_cxt[v_id] = x.value();
			}
			else
				PL_ERR("v_id {} has no left context", v_id);
		}

		for (pt::u32 j{start + len - 1}; j >= (start); j--) {
			pt::u32 v_id = rw->v_ids[j];

			if (pv_cmp::contains(right_cxt, v_id))
				continue;

			if (pv_cmp::contains(w_to_it, v_id))
				right_cxt[v_id] = w_to_it.at(v_id);
			else if (j + 1 < (start + len)) {
				if (pv_cmp::contains(w_to_it,
						     rw->v_ids[j + 1])) {
					y = w_to_it.at(rw->v_ids[j + 1]);
				}

				right_cxt[v_id] = y.value();
			}
			else
				PL_ERR("v_id {} has no right context", v_id);
		}
	};

	// fill sw with the first lap
	auto init = [&](const pt::slice &lap, const liteseq::ref_walk *rw)
	{
		auto [start, len] = lap.data();
		for (pt::u32 j{start}; j < (start + len); j++) {
			pt::u32 v_id = rw->v_ids[j];

			if (pv_cmp::contains(w_to_it, v_id))
				continue;

			auto it = sw.insert(sw.end(), v_id);
			w_to_it[v_id] = it;
		}
	};

	for (pt::u32 h_idx{}; h_idx < HAP_COUNT; h_idx++) {
		laps.clear();
		find_laps(g, h_idx, start_id, stop_id, laps);

		// if (dbg && !laps.empty())
		//	std::cerr << "->" << "\t";

		const liteseq::ref_walk *rw = g.get_ref_vec(h_idx)->walk;

		for (const auto &lap : laps) {
			auto [start, len] = lap.data();

			if (dbg) {
				std::cerr << "Lap: " << h_idx << " ";
				for (pt::u32 j{start}; j < (start + len); j++) {
					pt::u32 v_id = rw->v_ids[j];
					std::cerr << v_id << ",";
				}
				std::cerr << "\n";
			}

			pt::op_t<
				std::map<pt::u32, std::list<pt::u32>::iterator>>
				sw_cxt;

			if (sw.empty()) {
				init(lap, rw);
				if (dbg) {
					std::cerr << "~>\t"
						  << pu::concat_with(sw, ',')
						  << "\n";
				}
				continue;
			}
			else {
				comp_steps_cxt(lap, rw, sw_cxt);
			}

			const auto &[left_cxt, right_cxt] = sw_cxt;

			for (pt::u32 j{start}; j < (start + len); j++) {
				pt::u32 v_id = rw->v_ids[j];

				// if (dbg)
				//	std::cerr << v_id << ",";

				if (pv_cmp::contains(w_to_it, v_id))
					continue;

				/* v_id not in sw */

				// std::cerr << "trying " << v_id
				//	  << "left: " << *l_it
				//	  << " right: " << *r_it << "\n";

				auto l_it = left_cxt.at(v_id);
				auto r_it = right_cxt.at(v_id);

				// std::cerr << "trying " << v_id << " [" <<
				// *l_it
				//	  << ", " << *r_it << "]\n";

				auto insert_point =
					(std::distance(sw.begin(), l_it) <
					 std::distance(sw.begin(), r_it))
						? r_it
						: l_it;
				if (dbg) {
					std::cerr << "trying " << v_id << " ["
						  << *l_it << ", " << *r_it
						  << "] ip " << *insert_point
						  << "\n";
				}

				auto it = sw.insert(insert_point, v_id);
				w_to_it[v_id] = it;
			}

			if (dbg)
				std::cerr << "~>\t" << pu::concat_with(sw, ',')
					  << "\n";
		}
	}

	if (dbg) {
		INFO("sw: {}", rov.as_str());
		std::cerr << pu::concat_with(sw, ',') << "\n";
	}

	return sw;
}

std::list<pt::u32> gen_sort_old(const bd::VG &g, ir::RoV &rov,
				const pt::u32 HAP_COUNT, bool dbg)
{
	auto [l, r, route] = *rov.get_pvst_vtx()->get_route_params();
	auto [start_id, _] = l;
	auto [stop_id, __] = r;

	std::list<pt::u32> sw;
	std::map<pt::u32, std::list<pt::u32>::iterator> w_to_it;

	std::vector<pt::slice> laps;
	laps.reserve(1024);

	for (pt::u32 h_idx{}; h_idx < HAP_COUNT; h_idx++) {
		laps.clear();
		find_laps(g, h_idx, start_id, stop_id, laps);

		if (dbg && !laps.empty())
			std::cerr << "->" << "\t";

		const liteseq::ref_walk *rw = g.get_ref_vec(h_idx)->walk;

		for (const auto &lap : laps) {
			auto [start, len] = lap.data();

			for (pt::u32 j{start}; j < (start + len); j++) {
				pt::u32 v_id = rw->v_ids[j];

				if (dbg)
					std::cerr << v_id << ",";

				if (pv_cmp::contains(w_to_it, v_id)) {
					// move v_id to the back of the list
					auto it = w_to_it.at(v_id);
					sw.erase(it);
					it = sw.insert(sw.end(), v_id);
					w_to_it[v_id] = it;
				}
				else {
					auto it = sw.insert(sw.end(), v_id);
					w_to_it[v_id] = it;
				}
			}

			if (dbg) {
				std::cerr << "\n";
				std::cerr << "sw:\t" << pu::concat_with(sw, ',')
					  << "\n";
			}
		}
	}

	return sw;
}

pt::status_t enum_walks(const bd::VG &g, pvst::route_e route, idx_or_t src,
			idx_or_t snk, std::vector<ir::enhanced_walk> &walks,
			const std::string_view &rov_label)
{
	// default is source to sink
	dir_e ve_dir = OUT;
	dir_e nbr_dir = IN;
	idx_or_t start = src;
	idx_or_t end = snk;

	if (__builtin_expect(route == pvst::route_e::e2s, 0)) {
		ve_dir = IN;
		nbr_dir = OUT;
		start = snk;
		end = src;
	}

	std::queue<std::pair<pgt::walk_t, std::vector<pt::op_t<pt::u32>>>> q;
	q.push({{start}, {}});

	auto in_walk = [&](const pgt::walk_t &w, idx_or_t v) -> int
	{
		for (pt::u32 i = 0; i < w.size(); i++)
			if (w[i].v_id == v.v_id)
				return i;

		return -1;
	};

	auto to_id_walk = [&](const pgt::walk_t &w) -> pgt::walk_t
	{
		pgt::walk_t id_w;
		id_w.reserve(w.size());

		for (auto [v_idx, o] : w)
			id_w.push_back({g.v_idx_to_id(v_idx), o});

		return id_w;
	};

	while (!q.empty()) {
		// get the incoming vertices based on orientation
		auto [curr_w, cycles] = q.front();
		q.pop();

		auto w_len = static_cast<pt::u32>(curr_w.size());

		if (w_len > MAX_FLUBBLE_STEPS) {
			WARN("max steps reached for {}", rov_label);
			std::cerr << pgt::to_string(to_id_walk(curr_w)) << "\n";
			continue;
		}

		// discard a few walks
		if (q.size() > MAX_UNBLOCK_CTR) {
			continue;
		}

		const idx_or_t &curr = curr_w.back();
		if (curr == end) {
			// std::cerr << pgt::to_string(to_id_walk(curr_w)) <<
			// "\n";
			// // print cycles
			// for (auto [start_idx, end_idx] : cycles) {
			//	std::cerr << fmt::format("  cycle: {} -> {}\n",
			//				 start_idx, end_idx);
			// }

			walks.push_back({to_id_walk(curr_w), cycles});
			continue;
		}

		auto [v_idx, o] = curr;

		pgt::v_end_e ve = get_v_end(o, ve_dir);
		const bd::Vertex &v = g.get_vertex_by_idx(v_idx);
		const std::set<pt::idx_t> &nbr_edges = edges_at_end(v, ve);

		for (pt::u32 e_idx : nbr_edges) {
			const bd::Edge &e = g.get_edge(e_idx);
			auto [side, alt_idx] = e.get_other_vtx(v_idx, ve);
			idx_or_t nbr{alt_idx, get_or(side, nbr_dir)};

			int cycle_idx = in_walk(curr_w, nbr);
			if (cycle_idx >= 0) {
				// report cycle
				// std::cerr << "cycle detected: "
				//	  << curr_w.size() - 1 << " "
				//	  << cycle_idx << "\n";
				cycles.emplace_back(cycle_idx, w_len - 1);
				continue;
			}

			pgt::walk_t new_w = curr_w;
			new_w.emplace_back(nbr);

			// w_set.insert(nbr.v_id);

			q.emplace(new_w, cycles);
		}
	}

	return 0;
}

std::vector<pt::id_t> bfs_sort(const bd::VG &g, ir::RoV &rov)
{
	// Assume route parameters are already set.
	// Use structured bindings to unpack the pvst::route_params_t
	// object.
	auto [l, r, route] = *rov.get_pvst_vtx()->get_route_params();
	auto [start_id, start_o] = l;
	auto [stop_id, stop_o] = r;

	idx_or_t src = {g.v_id_to_idx(start_id), start_o};
	idx_or_t snk = {g.v_id_to_idx(stop_id), stop_o};

	ita::bfs::BfsTree t = comp_bfs_tree(g, route, src, snk);

	t.comp_depth();
	pt::status_t bfs_sort_res = t.sort();
	if (bfs_sort_res != 0)
		return {};

	return t.get_sorted();
}

pt::status_t find_walks(const bd::VG &g, ir::RoV &rov)
{
	bool dbg = rov.as_str() == ">1>4" ? true : false;

	// print each lap for each haplotype
	if (dbg) {
		auto [l, r, route] = *rov.get_pvst_vtx()->get_route_params();
		auto [start_id, _] = l;
		auto [stop_id, __] = r;

		std::vector<pt::slice> laps;
		laps.reserve(1024);

		for (pt::u32 h_idx{}; h_idx < g.get_hap_count(); h_idx++) {
			laps.clear();
			find_laps(g, h_idx, start_id, stop_id, laps);

			if (dbg && !laps.empty())
				std::cerr << "\n" << h_idx << "\t";

			const liteseq::ref_walk *rw =
				g.get_ref_vec(h_idx)->walk;

			for (const auto &lap : laps) {
				auto [start, len] = lap.data();

				for (pt::u32 j{start}; j < (start + len); j++) {
					pt::u32 v_id = rw->v_ids[j];

					if (dbg)
						std::cerr << v_id << ",";
				}
			}
		}

		std::cerr << "\n---------------------\n";
	}

	// Attempt to sort with BFS
	auto sorted_v_ids = bfs_sort(g, rov);
	if (!sorted_v_ids.empty()) {
		if (dbg)
			INFO("called 1");
		rov.add_sort_data(sorted_v_ids.begin(), sorted_v_ids.end());
		return 0; // Success
	}

	// Fallback: Generate sort data
	std::list<pt::u32> sw = gen_sort(g, rov, g.get_hap_count(), dbg);
	if (!sw.empty()) {
		if (dbg)
			INFO("called 2");
		rov.add_sort_data(sw.begin(), sw.end());
		return 0; // Success
	}

	return 1; // Failure
}
} // namespace povu::genomics::graph
