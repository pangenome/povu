#include "./smothered.hpp"
#include <string>


namespace povu::smothered {

struct fl_sls {
  pt::idx_t cn_v_idx;
  std::vector<pvst::Smothered> g_adj;
  std::vector<pvst::Smothered> s_adj;

  pt::idx_t size() const { return g_adj.size() + s_adj.size(); }

  // Constructor for fl_sls
  fl_sls(pt::idx_t cn_v_idx_)
    : cn_v_idx(cn_v_idx_),
      g_adj(std::vector<pvst::Smothered>{}),
      s_adj(std::vector<pvst::Smothered>{}) {}
};

pvst::bounds_t compute_bounds(const pst::Tree &st, pt::idx_t cn_st_idx, pt::idx_t sm_st_idx) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  if (st.is_desc(cn_st_idx, sm_st_idx)) {
    return pvst::bounds_t{cn_st_idx, sm_st_idx};
  }
  else {
    return pvst::bounds_t{sm_st_idx, cn_st_idx};
  }
}

bool is_nesting(const pst::Tree &st, const pvst::bounds_t &outer, const pvst::bounds_t &inner) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  return st.is_desc(outer.upper, inner.upper) && st.is_desc(inner.lower, outer.lower);
}


namespace g {

void trunk(const pst::Tree &st, const pvtr::Tree &pvst, const pvst::Concealed &ft_v,
           pt::idx_t cn_pvst_v_idx, const ptu::tree_meta &tm, std::vector<pvst::Smothered> &res) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  const std::vector<pt::idx_t> &depth = tm.depth;

  const pvst::Concealed &cn_v = ft_v;
  pt::idx_t fl_v_idx = cn_v.get_fl_idx();
  const pvst::Flubble fl_v = static_cast<const pvst::Flubble &>(pvst.get_vertex(fl_v_idx));
  pt::idx_t ai_st_idx = fl_v.get_ai();
  pt::idx_t zi = fl_v.get_zi();

  pt::idx_t sl_st_idx = ft_v.get_sl_st_idx();

  auto comp_e = [&](pt::idx_t st_v_idx) -> pgt::id_or_t {
    if (st.get_vertex(st_v_idx).type() == pst::v_type_e::l) {
      return { st.get_vertex(st_v_idx).g_v_id(), pgt::or_e::reverse };
    }
    else {
      return { st.get_vertex(st_v_idx).g_v_id(), pgt::or_e::forward };
    }
  };

  std::vector<pt::idx_t> mb;
  for (auto c_v_idx : st.get_children(sl_st_idx)) {
    std::unordered_set<pt::idx_t> c_br_srcs;
    pt::idx_t be_idx_ = pc::INVALID_IDX;
    for (auto be_idx : tm.get_brackets(c_v_idx)) {
      const pst::BackEdge &be = st.get_be(be_idx);
      pt::idx_t src = be.get_src();
      be_idx_ = be_idx;
      c_br_srcs.insert(src);
    }

    if (be_idx_ == pc::INVALID_IDX) {
      continue; // no brackets for this child
    }

    pt::idx_t tgt = st.get_be(be_idx_).get_tgt();
    pt::idx_t src = st.get_be(be_idx_).get_src();

    // compute bounds
    std::vector<pt::idx_t> p{zi, src};
    pt::idx_t brch_vtx = ptu::find_lca(tm, p);
    pvst::bounds_t bounds = pvst::bounds_t{brch_vtx, src};

    if (be_idx_ != pc::INVALID_IDX && c_br_srcs.size() == 1 && depth[tgt] > depth[ai_st_idx]) {
      pgt::id_or_t g = ft_v.get_cn_b();
      if (tm.get_brackets(src).empty()) {
        pgt::id_or_t e = comp_e(tgt);
        res.push_back(pvst::Smothered(g, e, cn_pvst_v_idx, false, tgt, pvst::cn_type_e::g, bounds));
      }
      else {
        pgt::id_or_t e = comp_e(src);
        res.push_back(pvst::Smothered(g, e, cn_pvst_v_idx, true, src, pvst::cn_type_e::g, bounds));
      }
    }
  }
  return;
}

void branch(const pst::Tree &st, const pvtr::Tree &pvst,
            const pvst::Concealed &ft_v, pt::idx_t cn_pvst_v_idx,
            const ptu::tree_meta &tm,
            std::vector<pvst::Smothered> &res) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  const pvst::Concealed &cn_v = ft_v;
  pt::idx_t fl_v_idx = cn_v.get_fl_idx();
  const pvst::Flubble fl_v = static_cast<const pvst::Flubble &>(pvst.get_vertex(fl_v_idx));
  pt::idx_t ai_st_idx = fl_v.get_ai();
  //pt::idx_t ai_st_idx = ft_v.get_ai();
  pt::idx_t sl_st_idx = ft_v.get_sl_st_idx();
  std::set<pt::idx_t> tgts = st.get_obe_tgt_v_idxs(sl_st_idx);

  if (tgts.size() != 1) {
    return;
  }

  pt::idx_t tgt = *tgts.begin();

  if (tgt != ai_st_idx) {
    return;
  }

  for (pt::idx_t be_idx : tm.get_brackets(sl_st_idx)) {
    const pst::BackEdge &be = st.get_be(be_idx);
    if(ai_st_idx == be.get_tgt()) {

      // there exists at least one bracket of the src that goes into sl_st_idx
      // TODO: what is be_idx_ here?
      for (auto be_idx_ : tm.get_brackets(be.get_src())){
        const pst::BackEdge &be = st.get_be(be_idx);
        if (be.get_tgt() == sl_st_idx) {
          pt::idx_t src = be.get_src();
          pgt::id_or_t e = (st.get_vertex(src).type() == pst::v_type_e::l)
                               ? pgt::id_or_t{st.get_vertex(src).g_v_id(), pgt::or_e::reverse}
                               : pgt::id_or_t{st.get_vertex(src).g_v_id(), pgt::or_e::forward};
          //pgt::id_or_t g = ft_v.get_end();
          pgt::id_or_t g = ft_v.get_cn_b();
          pvst::bounds_t bounds = compute_bounds(st, sl_st_idx, src);

          res.push_back(pvst::Smothered(g, e, cn_pvst_v_idx, true, src, pvst::cn_type_e::g, bounds));
          continue;
        }
      }

    }
  }

}
} // namespace g


namespace s {

pgt::id_or_t comp_w(const pst::Tree &st, pt::idx_t st_v_idx) {
  if (st.get_vertex(st_v_idx).type() == pst::v_type_e::l) {
    return {st.get_vertex(st_v_idx).g_v_id(), pgt::or_e::reverse};
  } else {
    return {st.get_vertex(st_v_idx).g_v_id(), pgt::or_e::forward};
  }
};

void trunk(const pst::Tree &st, const pvtr::Tree &pvst,
           const pvst::Concealed &ft_v, pt::idx_t cn_pvst_v_idx,
           const ptu::tree_meta &tm,
           std::vector<pvst::Smothered> &res) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  const std::vector<pt::idx_t> &depth = tm.depth;

  const pvst::Concealed &cn_v = ft_v;
  pt::idx_t fl_v_idx = cn_v.get_fl_idx();
  const pvst::Flubble fl_v = static_cast<const pvst::Flubble &>(pvst.get_vertex(fl_v_idx));
  pt::idx_t zi_st_idx = fl_v.get_zi();
  pt::idx_t sl_st_idx = ft_v.get_sl_st_idx();

  auto comp_w = [&](pt::idx_t st_v_idx) -> pgt::id_or_t {
    if (st.get_vertex(st_v_idx).type() == pst::v_type_e::l) {
      return { st.get_vertex(st_v_idx).g_v_id(), pgt::or_e::reverse };
    }
    else {
      return {st.get_vertex(st_v_idx).g_v_id(), pgt::or_e::forward};
    }
  };

  // key is LCA value is all the srcs
  std::map<pt::idx_t, std::vector<pt::idx_t>> lca_map;

  std::set<pt::idx_t> srcs = st.get_ibe_src_v_idxs(sl_st_idx);
  for (auto src: srcs) {
    if (depth[src] >= depth[zi_st_idx]) {
      continue;
    }

    std::vector<pt::idx_t> p {zi_st_idx, src};
    pt::idx_t lca = ptu::find_lca(tm, p);
    if (lca == src) {
      // also lca <= zi_st_idx
      continue; // lca is in the trunk
    }
    lca_map[lca].push_back(src);
  }

  for (auto [lca, srcs] : lca_map) {
    if (srcs.size() == 1) {
      pt::idx_t src = srcs[0];
      pgt::id_or_t w = comp_w(src);
      pgt::id_or_t s = ft_v.get_cn_b();
      pvst::bounds_t bounds = {lca, src};

      res.push_back(pvst::Smothered(s, w, cn_pvst_v_idx, false, src, pvst::cn_type_e::s, bounds));
    }
  }
}

void branch(const pst::Tree &st, const pvtr::Tree &pvst,
            const pvst::Concealed &ft_v, pt::idx_t cn_pvst_v_idx,
            const ptu::tree_meta &tm,
            std::vector<pvst::Smothered> &res) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  const pvst::Concealed &cn_v = ft_v;
  pt::idx_t fl_v_idx = cn_v.get_fl_idx();
  const pvst::Flubble fl_v = static_cast<const pvst::Flubble &>(pvst.get_vertex(fl_v_idx));
  pt::idx_t zi_st_idx = fl_v.get_zi();
  pt::idx_t sl_st_idx = ft_v.get_sl_st_idx();

  std::set<pt::idx_t> srcs = st.get_ibe_src_v_idxs(sl_st_idx);

  if (srcs.empty()) {
    return; // no srcs
  }

  for (pt::idx_t src_v_idx : srcs) {
    for (auto src_ : st.get_ibe_src_v_idxs(src_v_idx)) {
      pgt::id_or_t w = comp_w(st, src_);
      pgt::id_or_t s = ft_v.get_cn_b();
      pvst::bounds_t bounds = compute_bounds(st, src_, sl_st_idx);

      res.push_back(pvst::Smothered(s, w, cn_pvst_v_idx, true, src_, pvst::cn_type_e::s, bounds));
    }

    for (auto tgt_ : st.get_obe_tgt_v_idxs(src_v_idx)) {
      if (tgt_ == zi_st_idx) {
        continue; // skip the trunk
      }

      pgt::id_or_t w = comp_w(st, tgt_);
      pgt::id_or_t s = ft_v.get_cn_b();
      pvst::bounds_t bounds = compute_bounds(st, tgt_, sl_st_idx);

      res.push_back(pvst::Smothered(s, w, cn_pvst_v_idx, false, tgt_, pvst::cn_type_e::s, bounds));
    }
  }
}

} // namespace s

void nest(const pst::Tree &st, pvtr::Tree &pvst, const ptu::tree_meta &tm,
          pt::idx_t cn_pvst_v_idx, pt::idx_t smo_pvst_v_idx) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  const pvst::Smothered &smo_v = static_cast<const pvst::Smothered &>(pvst.get_vertex(smo_pvst_v_idx));

  std::vector<pt::idx_t> ch = pvst.get_children(cn_pvst_v_idx);
  for (pt::idx_t c_v_idx : ch) {
    pvst::VertexBase &pvst_v = pvst.get_vertex_mut(c_v_idx);
    pvst::bounds_t bounds = {pc::INVALID_IDX, pc::INVALID_IDX};

    std::cerr << fn_name << " " << pvst_v.as_str() << "\n";

    switch (pvst_v.get_type()) {
    case pvst::vt_e::flubble:
    case pvst::vt_e::tiny:
    case pvst::vt_e::parallel:
      bounds = static_cast<const pvst::Flubble &>(pvst_v).get_bounds();
      break;
    case pvst::vt_e::slubble: // concealed
      bounds = static_cast<const pvst::Concealed &>(pvst_v).get_bounds();
      break;
    default:
      break;
    }

    if (bounds.upper == pc::INVALID_IDX || bounds.lower == pc::INVALID_IDX) {
      continue;
    }


    if (is_nesting(st, smo_v.get_bounds(), bounds)) {
      pvst.del_edge(cn_pvst_v_idx, c_v_idx);
      pvst.add_edge(smo_pvst_v_idx, c_v_idx);
    }
  }
}

void add_smothered(const pst::Tree &st, pvtr::Tree &pvst,
                   const ptu::tree_meta &tm, const std::vector<fl_sls> &al_smo) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  for (const fl_sls &smo : al_smo) {

     for (const pvst::Smothered &g_adj_v : smo.g_adj) {
      std::cerr << fn_name << " adding smothered [e,g] vertex: " << g_adj_v.as_str() << "\n";
      // use a ref and move?

      pt::idx_t g_adj_v_idx = pvst.add_vertex(g_adj_v);
      pvst.add_edge(smo.cn_v_idx, g_adj_v_idx);
      nest(st, pvst, tm, smo.cn_v_idx, g_adj_v_idx);
    }

    for (const pvst::Smothered &s_adj_v : smo.s_adj) {
      std::cerr << fn_name << " adding smothered [s,w] vertex: " << s_adj_v.as_str() << "\n";
      // use a ref and move?
      pt::idx_t s_adj_v_idx = pvst.add_vertex(s_adj_v);
      pvst.add_edge(smo.cn_v_idx, s_adj_v_idx);
      nest(st, pvst, tm, smo.cn_v_idx, s_adj_v_idx);
    }
  }

  return;
}

void find_smothered(const pst::Tree &st, pvtr::Tree &ft, const ptu::tree_meta &tm) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  pvtr::Tree &pvst = ft;

  std::vector<fl_sls> all_smo;

  for (pt::idx_t ft_v_idx{}; ft_v_idx < ft.vtx_count(); ft_v_idx++) {
     const pvst::VertexBase &pvst_v = ft.get_vertex(ft_v_idx);

     if (pvst_v.get_type() != pvst::vt_e::slubble) {
       continue;
     }

     const pvst::Concealed &cn_v = static_cast<const pvst::Concealed &>(pvst_v);

     fl_sls smo {ft_v_idx};
     switch (cn_v.get_sl_type()) {
     case pvst::sl_type_e::ai_trunk:
       g::trunk(st, pvst, cn_v, ft_v_idx, tm, smo.g_adj);
       break;
     case pvst::sl_type_e::ai_branch:
       g::branch(st, pvst, cn_v, ft_v_idx, tm, smo.g_adj);
       break;
     case pvst::sl_type_e::zi_trunk:
       s::trunk(st, pvst, cn_v, ft_v_idx, tm, smo.s_adj);
       break;
     case pvst::sl_type_e::zi_branch:
       s::branch(st, pvst, cn_v, ft_v_idx, tm, smo.s_adj);
       break;
     default:
       std::cerr << fn_name << " unknown slubble type: " << ft_v_idx << "\n";
       continue; // skip this vertex
     }

    if (smo.size() > 0) {
      all_smo.push_back(smo);
    }
  }

  add_smothered(st, ft, tm, all_smo);

}

} // namespace povu::smothered
