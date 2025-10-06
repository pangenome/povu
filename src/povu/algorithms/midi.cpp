#include "povu/algorithms/midi.hpp"

#include <exception> // for exception
#include <map>	     // for map
#include <memory>    // for make_unique
#include <stdexcept> // for runtime_error
#include <string>    // for basic_string, string
#include <utility>   // for get, pair
#include <vector>    // for vector

#include "fmt/core.h"		  // for format
#include "povu/common/compat.hpp" // for format, pv_cmp
#include "povu/common/core.hpp"	  // for pt, idx_t
#include "povu/common/log.hpp"	  // for WARN, ERR

namespace povu::midi
{

// add midi bubbles to the PVST
void add_midi(const ptu::tree_meta &tm,
	      std::map<pt::idx_t, std::vector<pvst::MidiBubble>> &midis,
	      pvst::Tree &pvst)
{
	const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

	const std::vector<pt::idx_t> &depth = tm.depth;

	for (auto &[ft_v_idx, bubs] : midis) {
		std::vector<pt::idx_t> ch = pvst.get_children(ft_v_idx);

		for (const pvst::MidiBubble &b : bubs) {
			auto [md_upper, md_lower] =
				b.get_bounds(); // get bounds of the midibubble

			// TODO: [A] investigate why this happens
			if (md_upper == pc::INVALID_IDX ||
			    md_lower == pc::INVALID_IDX) {
				continue;
			}

			pt::idx_t m_v_idx = pvst.add_vertex(b);
			pvst.add_edge(ft_v_idx, m_v_idx);

			// nest
			for (pt::idx_t c_v_idx : ch) {
				const pvst::VertexBase &c_v =
					pvst.get_vertex(c_v_idx);

				if (c_v.get_fam() != pvst::vt_e::flubble) {
					continue;
				}

				const pvst::Flubble &fl_v =
					static_cast<const pvst::Flubble &>(c_v);
				auto [fl_upper, fl_lower] = fl_v.get_bounds();
				if ((depth[md_upper] < depth[fl_upper]) &&
				    (depth[md_lower] > depth[fl_lower])) {
					// flubble is nested in the midibubble
					pvst.del_edge(ft_v_idx, c_v_idx);
					pvst.add_edge(m_v_idx, c_v_idx);
				}
			}
		}
	}
}

pvst::MidiBubble gen_midi_bub(const pvst::Tree &pvst,
			      const std::vector<pt::idx_t> &c_bubs)
{
	const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

	pt::idx_t fst = c_bubs[0]; // first
	pt::idx_t snd = c_bubs[1]; // second

	auto get_cn = [&](pt::idx_t cn_v_idx_pvst) -> pvst::Concealed
	{
		const pvst::VertexBase &pvst_v = pvst.get_vertex(cn_v_idx_pvst);
		const pvst::Concealed &cn_v =
			static_cast<const pvst::Concealed &>(pvst_v);
		return cn_v;
	};

	pgt::id_or_t g;
	pt::idx_t g_pvst_idx{pc::INVALID_IDX};
	pgt::id_or_t s;
	pt::idx_t s_pvst_idx{pc::INVALID_IDX};

	pvst::Concealed f_cn_v = get_cn(fst);
	if ((f_cn_v.get_sl_type() == pvst::cl_e::ai_branch) ||
	    (f_cn_v.get_sl_type() == pvst::cl_e::ai_trunk)) {
		// is g
		g = f_cn_v.get_cn_b();
		g_pvst_idx = fst;
	}
	else {
		s = f_cn_v.get_cn_b();
		s_pvst_idx = fst;
	}

	pvst::Concealed snd_cn_v = get_cn(snd);
	if ((snd_cn_v.get_sl_type() == pvst::cl_e::ai_branch) ||
	    (snd_cn_v.get_sl_type() == pvst::cl_e::ai_trunk)) {
		// is g
		g = snd_cn_v.get_cn_b();
		g_pvst_idx = snd;
	}
	else {
		s = snd_cn_v.get_cn_b();
		s_pvst_idx = snd;
	}

	// Validate indices before creating MidiBubble
	if (g_pvst_idx == pc::INVALID_IDX || s_pvst_idx == pc::INVALID_IDX) {
		ERR("{} Invalid MidiBubble indices: g_pvst_idx={}, "
		    "s_pvst_idx={}, fst={}, snd={}",
		    fn_name, g_pvst_idx, s_pvst_idx, fst, snd);
		throw std::runtime_error(
			"Cannot create MidiBubble with invalid indices");
	}

	// Ensure indices are within reasonable bounds
	pt::idx_t max_valid_idx = pvst.vtx_count();
	if (g_pvst_idx >= max_valid_idx || s_pvst_idx >= max_valid_idx) {
		ERR("{} MidiBubble indices out of bounds: g_pvst_idx={}, "
		    "s_pvst_idx={}, max_valid={}",
		    fn_name, g_pvst_idx, s_pvst_idx, max_valid_idx);
		throw std::runtime_error(
			"MidiBubble indices exceed PVST vertex count");
	}

	return pvst::MidiBubble::create(g_pvst_idx, g, s_pvst_idx, s,
					pvst::route_e::s2e);
}

// find midibubbles in the flubble
std::vector<pvst::MidiBubble> handle_fl(const pst::Tree &st,
					const pvst::Tree &pvst,
					const pvst::Flubble &fl_v,
					const std::vector<pt::idx_t> &c_bubs)
{

	std::vector<pvst::MidiBubble> res;

	pt::idx_t zi = fl_v.get_zi();

	// is a descendant
	auto is_desc = [&](pt::idx_t p_v_idx, pt::idx_t c_v_idx) -> bool
	{
		return ((st.get_vertex(p_v_idx).pre_order() <
			 st.get_vertex(c_v_idx).pre_order()) &&
			(st.get_vertex(p_v_idx).post_order() >
			 st.get_vertex(c_v_idx).post_order()));
	};

	// is v_idx in the trunk? true if v_idx is an ancestor of zi, then it is
	// in the trunk
	auto in_trunk = [&](pt::idx_t v_idx) -> bool
	{
		return v_idx < zi;
	};

	std::vector<pt::idx_t> trunk_c_bubs;

	// key is child and value the set of concealed bubbles in the branch
	std::map<pt::idx_t, std::vector<pt::idx_t>> branch_c_bubs;

	std::vector<pt::idx_t> c; // children of z_i \ z_x
	for (pt::idx_t e_idx : st.get_child_edge_idxs(zi)) {
		const pst::Edge &e = st.get_tree_edge(e_idx);
		if (e.get_color() == pgt::color_e::black) {
			continue;
		}
		c.push_back(e.get_child_v_idx());
	}

	for (pt::idx_t cn_v_idx_pvst : c_bubs) {
		const pvst::VertexBase &pvst_v = pvst.get_vertex(cn_v_idx_pvst);
		const pvst::Concealed &cn_v =
			static_cast<const pvst::Concealed &>(pvst_v);
		pt::idx_t g_or_s_st_idx = cn_v.get_sl_st_idx();

		if (in_trunk(g_or_s_st_idx)) {
			trunk_c_bubs.push_back(cn_v_idx_pvst);
		}
		else {
			for (auto p_v_idx : c) {
				if (is_desc(p_v_idx, g_or_s_st_idx)) {
					branch_c_bubs[p_v_idx].push_back(
						g_or_s_st_idx);
				}
			}
		}
	}

	if (trunk_c_bubs.size() == 2) {
		// found midibubble in trunk
		// gen midibubble
		try {
			res.push_back(gen_midi_bub(pvst, trunk_c_bubs));
		}
		catch (const std::exception &e) {
			WARN("Failed to generate trunk MidiBubble: {}",
			     e.what());
		}
	}

	for (auto [c_idx, v] : branch_c_bubs) {
		if (v.size() == 2) {
			// found midibubble in branch
			// gen midibubble
			try {
				res.push_back(gen_midi_bub(pvst, v));
			}
			catch (const std::exception &e) {
				WARN("Failed to generate branch MidiBubble for "
				     "c_idx {}: {}",
				     c_idx, e.what());
			}
		}
	}

	return res;
}

void find_midi(const pst::Tree &st, pvst::Tree &pvst, const ptu::tree_meta &tm)
{
	const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

	std::map<pt::idx_t, std::vector<pvst::MidiBubble>> x;
	for (pt::idx_t ft_v_idx{}; ft_v_idx < pvst.vtx_count(); ft_v_idx++) {
		const pvst::VertexBase &pvst_v = pvst.get_vertex(ft_v_idx);

		if (pvst_v.get_fam() != pvst::vt_e::flubble) {
			continue;
		}

		std::vector<pt::idx_t> c_bubs; // the concealed bubbles
		for (auto c_v_idx_pvst : pvst.get_children(ft_v_idx)) {
			const pvst::VertexBase &c_v =
				pvst.get_vertex(c_v_idx_pvst);
			if (c_v.get_fam() == pvst::vt_e::concealed) {
				c_bubs.push_back(c_v_idx_pvst);
			}
		}

		if (c_bubs.size() < 2) {
			continue;
		}

		const pvst::Flubble &fl_v =
			static_cast<const pvst::Flubble &>(pvst_v);

		try {
			std::vector<pvst::MidiBubble> res =
				handle_fl(st, pvst, fl_v, c_bubs);
			if (!res.empty()) {
				x[ft_v_idx] = res;
			}
		}
		catch (const std::exception &e) {
			WARN("{} Failed to handle flubble {}: {}", fn_name,
			     ft_v_idx, e.what());
			continue;
		}
	}

	add_midi(tm, x, pvst);
}
} // namespace povu::midi
