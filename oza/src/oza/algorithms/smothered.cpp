#include "oza/algorithms/smothered.hpp"

#include <iostream>	 // for basic_ostream, operator<<, basic_ios
#include <map>		 // for map
#include <optional>	 // for optional, nullopt, nullopt_t
#include <set>		 // for set
#include <unordered_set> // for unordered_set
#include <utility>	 // for get, pair
#include <vector>	 // for vector

#include <log/location.hpp>	 // for LOG_HERE
#include <quilt/graph_types.hpp> // for v_end_e, side_n_id_t, side_n_idx_t
#include <quilt/shim.hpp>	 // for format
#include <quilt/types.hpp>	 // for qt

namespace oza::smothered
{
namespace pgt = quilt::types::graph;

struct fl_sls {
	qt::idx_t cn_v_idx;
	std::vector<pvst::Smothered> g_adj;
	std::vector<pvst::Smothered> s_adj;

	qt::idx_t size() const
	{
		return g_adj.size() + s_adj.size();
	}

	// Constructor for fl_sls
	fl_sls(qt::idx_t cn_v_idx_)
	    : cn_v_idx(cn_v_idx_), g_adj(std::vector<pvst::Smothered>{}),
	      s_adj(std::vector<pvst::Smothered>{})
	{}
};

pvst::bounds_t compute_bounds(const pst::Tree &st, qt::idx_t cn_st_idx,
			      qt::idx_t sm_st_idx)
{

	if (st.is_desc(cn_st_idx, sm_st_idx)) {
		return pvst::bounds_t{cn_st_idx, sm_st_idx};
	}
	else {
		return pvst::bounds_t{sm_st_idx, cn_st_idx};
	}
}

bool is_nesting(const pst::Tree &st, const pvst::bounds_t &outer,
		const pvst::bounds_t &inner)
{

	return st.is_desc(outer.upper, inner.upper) &&
	       st.is_desc(inner.lower, outer.lower);
}

namespace g
{

void trunk(const pst::Tree &st, const pvst::Tree &pvst,
	   const pvst::Concealed &ft_v, qt::idx_t cn_pvst_v_idx,
	   const ptu::tree_meta &tm, std::vector<pvst::Smothered> &res)
{

	const std::vector<qt::idx_t> &depth = tm.depth;

	const pvst::Concealed &cn_v = ft_v;
	qt::idx_t fl_v_idx = cn_v.get_fl_idx();
	const pvst::Flubble fl_v =
		static_cast<const pvst::Flubble &>(pvst.get_vertex(fl_v_idx));
	qt::idx_t ai_st_idx = fl_v.get_ai();
	qt::idx_t zi = fl_v.get_zi();

	qt::idx_t sl_st_idx = ft_v.get_sl_st_idx();

	auto comp_e = [&](qt::idx_t st_v_idx) -> pgt::id_or_t
	{
		if (st.get_vertex(st_v_idx).type() == pst::v_type_e::l) {
			return {st.get_vertex(st_v_idx).g_v_id(),
				pgt::or_e::reverse};
		}
		else {
			return {st.get_vertex(st_v_idx).g_v_id(),
				pgt::or_e::forward};
		}
	};

	std::vector<qt::idx_t> mb;
	for (auto c_v_idx : st.get_children(sl_st_idx)) {
		std::unordered_set<qt::idx_t> c_br_srcs;
		qt::idx_t be_idx_ = pc::INVALID_IDX;
		for (auto be_idx : tm.get_brackets(c_v_idx)) {
			const pst::BackEdge &be = st.get_be(be_idx);
			qt::idx_t src = be.get_src();
			be_idx_ = be_idx;
			c_br_srcs.insert(src);
		}

		if (be_idx_ == pc::INVALID_IDX) {
			continue; // no brackets for this child
		}

		qt::idx_t tgt = st.get_be(be_idx_).get_tgt();
		qt::idx_t src = st.get_be(be_idx_).get_src();

		// compute bounds
		std::vector<qt::idx_t> p{zi, src};
		qt::idx_t brch_vtx = ptu::find_lca(tm, p);
		pvst::bounds_t bounds = pvst::bounds_t{brch_vtx, src};

		if (be_idx_ != pc::INVALID_IDX && c_br_srcs.size() == 1 &&
		    depth[tgt] > depth[ai_st_idx]) {
			pgt::id_or_t g = ft_v.get_cn_b();
			if (tm.get_brackets(src).empty()) {
				pgt::id_or_t e = comp_e(tgt);
				res.push_back(pvst::Smothered::create(
					g, e, cn_pvst_v_idx, false, tgt,
					pvst::cb_e::g, bounds,
					pvst::rt_e::s2e));
			}
			else {
				pgt::id_or_t e = comp_e(src);
				res.push_back(pvst::Smothered::create(
					g, e, cn_pvst_v_idx, true, src,
					pvst::cb_e::g, bounds,
					pvst::rt_e::s2e));
			}
		}
	}
	return;
}

void branch(const pst::Tree &st, const pvst::Tree &pvst,
	    const pvst::Concealed &ft_v, qt::idx_t cn_pvst_v_idx,
	    const ptu::tree_meta &tm, std::vector<pvst::Smothered> &res)
{

	const pvst::Concealed &cn_v = ft_v;
	qt::idx_t fl_v_idx = cn_v.get_fl_idx();
	const pvst::Flubble fl_v =
		static_cast<const pvst::Flubble &>(pvst.get_vertex(fl_v_idx));
	qt::idx_t ai_st_idx = fl_v.get_ai();
	// qt::idx_t ai_st_idx = ft_v.get_ai();
	qt::idx_t sl_st_idx = ft_v.get_sl_st_idx();
	std::set<qt::idx_t> tgts = st.get_obe_tgt_v_idxs(sl_st_idx);

	if (tgts.size() != 1) {
		return;
	}

	qt::idx_t tgt = *tgts.begin();

	if (tgt != ai_st_idx) {
		return;
	}

	for (qt::idx_t be_idx : tm.get_brackets(sl_st_idx)) {
		const pst::BackEdge &be = st.get_be(be_idx);
		if (ai_st_idx == be.get_tgt()) {

			// there exists at least one bracket of the src that
			// goes into sl_st_idx
			// TODO: what is be_idx_ here?
			//  [A] fix broken loop here
			for (auto be_idx_ : tm.get_brackets(be.get_src())) {
				const pst::BackEdge &be_ = st.get_be(
					be_idx_); // should be get bracket?
				if (be_.get_tgt() == sl_st_idx) {
					qt::idx_t src = be.get_src();
					pgt::id_or_t e =
						(st.get_vertex(src).type() ==
						 pst::v_type_e::l)
							? pgt::id_or_t{st.get_vertex(
										 src)
									       .g_v_id(),
								       pgt::or_e::
									       reverse}
							: pgt::id_or_t{
								  st.get_vertex(
									    src)
									  .g_v_id(),
								  pgt::or_e::
									  forward};
					// pgt::id_or_t g = ft_v.get_end();
					pgt::id_or_t g = ft_v.get_cn_b();
					pvst::bounds_t bounds = compute_bounds(
						st, sl_st_idx, src);

					res.push_back(pvst::Smothered::create(
						g, e, cn_pvst_v_idx, true, src,
						pvst::cb_e::g, bounds,
						pvst::rt_e::s2e));
					continue;
				}
			}
		}
	}
}
} // namespace g

namespace s
{

pgt::id_or_t comp_w(const pst::Tree &st, qt::idx_t st_v_idx)
{
	if (st.get_vertex(st_v_idx).type() == pst::v_type_e::l) {
		return {st.get_vertex(st_v_idx).g_v_id(), pgt::or_e::reverse};
	}
	else {
		return {st.get_vertex(st_v_idx).g_v_id(), pgt::or_e::forward};
	}
};

void trunk(const pst::Tree &st, const pvst::Tree &pvst,
	   const pvst::Concealed &ft_v, qt::idx_t cn_pvst_v_idx,
	   const ptu::tree_meta &tm, std::vector<pvst::Smothered> &res)
{

	const std::vector<qt::idx_t> &depth = tm.depth;

	const pvst::Concealed &cn_v = ft_v;
	qt::idx_t fl_v_idx = cn_v.get_fl_idx();
	const pvst::Flubble fl_v =
		static_cast<const pvst::Flubble &>(pvst.get_vertex(fl_v_idx));
	qt::idx_t zi_st_idx = fl_v.get_zi();
	qt::idx_t sl_st_idx = ft_v.get_sl_st_idx();

	auto comp_w = [&](qt::idx_t st_v_idx) -> pgt::id_or_t
	{
		if (st.get_vertex(st_v_idx).type() == pst::v_type_e::l) {
			return {st.get_vertex(st_v_idx).g_v_id(),
				pgt::or_e::reverse};
		}
		else {
			return {st.get_vertex(st_v_idx).g_v_id(),
				pgt::or_e::forward};
		}
	};

	// key is LCA value is all the srcs
	std::map<qt::idx_t, std::vector<qt::idx_t>> lca_map;

	for (auto src : st.get_ibe_src_v_idxs(sl_st_idx)) {
		if (depth[src] >= depth[zi_st_idx]) {
			continue;
		}

		std::vector<qt::idx_t> p{zi_st_idx, src};
		qt::idx_t lca = ptu::find_lca(tm, p);
		if (lca == src) {
			// also lca <= zi_st_idx
			continue; // lca is in the trunk
		}
		lca_map[lca].push_back(src);
	}

	for (auto [lca, srcs] : lca_map) {
		if (srcs.size() == 1) {
			qt::idx_t src = srcs[0];
			pgt::id_or_t w = comp_w(src);
			pgt::id_or_t s = ft_v.get_cn_b();
			pvst::bounds_t bounds = {lca, src};

			res.push_back(pvst::Smothered::create(
				s, w, cn_pvst_v_idx, false, src, pvst::cb_e::s,
				bounds, pvst::rt_e::e2s));
		}
	}
}

void branch(const pst::Tree &st, const pvst::Tree &pvst,
	    const pvst::Concealed &ft_v, qt::idx_t cn_pvst_v_idx,
	    std::vector<pvst::Smothered> &res)
{

	const pvst::Concealed &cn_v = ft_v;
	qt::idx_t fl_v_idx = cn_v.get_fl_idx();
	const pvst::Flubble fl_v =
		static_cast<const pvst::Flubble &>(pvst.get_vertex(fl_v_idx));
	qt::idx_t zi_st_idx = fl_v.get_zi();
	qt::idx_t sl_st_idx = ft_v.get_sl_st_idx();

	std::set<qt::idx_t> srcs = st.get_ibe_src_v_idxs(sl_st_idx);

	if (srcs.empty()) {
		return; // no srcs
	}

	for (qt::idx_t src_v_idx : srcs) {
		for (auto src_ : st.get_ibe_src_v_idxs(src_v_idx)) {
			pgt::id_or_t w = comp_w(st, src_);
			pgt::id_or_t s = ft_v.get_cn_b();
			pvst::bounds_t bounds =
				compute_bounds(st, src_, sl_st_idx);

			res.push_back(pvst::Smothered::create(
				s, w, cn_pvst_v_idx, true, src_, pvst::cb_e::s,
				bounds, pvst::rt_e::e2s));
		}

		for (auto tgt_ : st.get_obe_tgt_v_idxs(src_v_idx)) {
			if (tgt_ == zi_st_idx) {
				continue; // skip the trunk
			}

			pgt::id_or_t w = comp_w(st, tgt_);
			pgt::id_or_t s = ft_v.get_cn_b();
			pvst::bounds_t bounds =
				compute_bounds(st, tgt_, sl_st_idx);

			res.push_back(pvst::Smothered::create(
				s, w, cn_pvst_v_idx, false, tgt_, pvst::cb_e::s,
				bounds, pvst::rt_e::e2s));
		}
	}
}

} // namespace s

[[nodiscard]] inline std::optional<pvst::bounds_t>
maybe_bounds(const pvst::VertexBase &v) noexcept
{
	if (auto f = dynamic_cast<const pvst::Flubble *>(&v))
		return f->get_bounds();
	if (auto c = dynamic_cast<const pvst::Concealed *>(&v))
		return c->get_bounds();
	return std::nullopt;
}

void nest(const pst::Tree &st, pvst::Tree &pvst, qt::idx_t cn_pvst_v_idx,
	  qt::idx_t smo_pvst_v_idx)
{

	const pvst::Smothered &smo_v = static_cast<const pvst::Smothered &>(
		pvst.get_vertex(smo_pvst_v_idx));

	for (qt::idx_t c_v_idx : pvst.get_children(cn_pvst_v_idx)) {

		const pvst::VertexBase &pvst_v = pvst.get_vertex(c_v_idx);

		if (auto b = maybe_bounds(pvst_v);
		    b && is_nesting(st, smo_v.get_bounds(), *b)) {
			pvst.del_edge(cn_pvst_v_idx, c_v_idx);
			pvst.add_edge(smo_pvst_v_idx, c_v_idx);
		}
	}
}

void add_smothered(const pst::Tree &st, pvst::Tree &pvst,
		   const std::vector<fl_sls> &al_smo)
{

	for (const fl_sls &smo : al_smo) {

		for (const pvst::Smothered &g_adj_v : smo.g_adj) {
			// std::cerr << fn_name << " adding smothered [e,g]
			// vertex: " << g_adj_v.as_str() << "\n";
			// use a ref and move?

			qt::idx_t g_adj_v_idx = pvst.add_vertex(g_adj_v);
			pvst.add_edge(smo.cn_v_idx, g_adj_v_idx);
			nest(st, pvst, smo.cn_v_idx, g_adj_v_idx);
		}

		for (const pvst::Smothered &s_adj_v : smo.s_adj) {
			// std::cerr << fn_name << " adding smothered [s,w]
			// vertex: " << s_adj_v.as_str() << "\n";
			//  use a ref and move?
			qt::idx_t s_adj_v_idx = pvst.add_vertex(s_adj_v);
			pvst.add_edge(smo.cn_v_idx, s_adj_v_idx);
			nest(st, pvst, smo.cn_v_idx, s_adj_v_idx);
		}
	}

	return;
}

void find_smothered(const pst::Tree &st, pvst::Tree &ft,
		    const ptu::tree_meta &tm)
{

	pvst::Tree &pvst = ft;

	std::vector<fl_sls> all_smo;

	for (qt::idx_t ft_v_idx{}; ft_v_idx < ft.vtx_count(); ft_v_idx++) {
		const pvst::VertexBase &pvst_v = ft.get_vertex(ft_v_idx);

		if (pvst_v.get_fam() != pvst::vf_e::concealed) {
			continue;
		}

		const pvst::Concealed &cn_v =
			static_cast<const pvst::Concealed &>(pvst_v);

		fl_sls smo{ft_v_idx};
		switch (cn_v.get_sl_type()) {
		case pvst::cl_e::ai_trunk:
			g::trunk(st, pvst, cn_v, ft_v_idx, tm, smo.g_adj);
			break;
		case pvst::cl_e::ai_branch:
			g::branch(st, pvst, cn_v, ft_v_idx, tm, smo.g_adj);
			break;
		case pvst::cl_e::zi_trunk:
			s::trunk(st, pvst, cn_v, ft_v_idx, tm, smo.s_adj);
			break;
		case pvst::cl_e::zi_branch:
			s::branch(st, pvst, cn_v, ft_v_idx, smo.s_adj);
			break;
		default:
			std::cerr << LOG_HERE
				  << " unknown slubble type: " << ft_v_idx
				  << "\n";
			continue; // skip this vertex
		}

		if (smo.size() > 0) {
			all_smo.push_back(smo);
		}
	}

	add_smothered(st, ft, all_smo);
}

} // namespace oza::smothered
