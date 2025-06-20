#include "./smothered.hpp"


namespace povu::smothered {

struct fl_sls {
  pt::idx_t fl_v_idx;
  std::vector<pvst::Vertex> g_adj;
  std::vector<pvst::Vertex> s_adj;

  pt::idx_t size() const { return g_adj.size() + s_adj.size(); }

  fl_sls(pt::idx_t fl_v_idx)
    : g_adj(std::vector<pvst::Vertex>{}),
      s_adj(std::vector<pvst::Vertex>{}) {}
  };

namespace g {

void trunk(const pst::Tree &st, const pvst::Vertex &ft_v,
           const ptu::tree_meta &tm, std::vector<pvst::Vertex> &res) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  //std::vector<pvst::Vertex> res;
  //std::vector<pt::idx_t> x;

  const std::vector<pt::idx_t> &depth = tm.depth;

  pt::idx_t ai_st_idx = ft_v.get_ai();
  //pt::idx_t zi = ft_v.get_zi();
  pt::idx_t sl_st_idx = ft_v.get_sl_st_idx();

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

    // std::cerr << "ai: " << ai_st_idx << " tgt " << tgt << "\n";

    if (be_idx_ != pc::INVALID_IDX && c_br_srcs.size() == 1 && depth[tgt] > depth[ai_st_idx]) {

      //pgt::id_or_t e;
      if (tm.get_brackets(src).empty()) {
        pgt::id_or_t e =
            (st.get_vertex(tgt).type() == pst::v_type_e::l)
                ? pgt::id_or_t{st.get_vertex(tgt).g_v_id(), pgt::or_e::reverse}
                : pgt::id_or_t{st.get_vertex(tgt).g_v_id(), pgt::or_e::forward};
        pgt::id_or_t g = ft_v.get_end();
        res.push_back(pvst::Vertex::make_smothered(e, g, tgt, ft_v));
      }
      else {
        pgt::id_or_t e =
            (st.get_vertex(src).type() == pst::v_type_e::l)
                ? pgt::id_or_t{st.get_vertex(src).g_v_id(), pgt::or_e::reverse}
                : pgt::id_or_t{st.get_vertex(src).g_v_id(), pgt::or_e::forward};
        pgt::id_or_t g = ft_v.get_end();
        res.push_back(pvst::Vertex::make_smothered(e, g, src,ft_v));
      }

    }
  }

  // for (const pvst::Vertex &v : res) {
  //   std::cerr << fn_name << " found [1] (a): " << v.as_str() << "\n";

  //   for (auto be_idx : tm.get_brackets(sl_st_idx)) {
  //     const pst::BackEdge &be = st.get_be(be_idx);

  //     if (be.get_tgt() != ai_st_idx) {
  //       continue; // not a back edge to the ai_st_idx
  //     }
  //     //pt::idx_t src = be.get_src();
  //   }
  //}
}

void branch(const pst::Tree &st, const pvst::Vertex &ft_v,
            const ptu::tree_meta &tm, std::vector<pvst::Vertex> &res) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  //std::vector<pt::idx_t> e;
  //g_branches(st, ft_v, tm, e);

  const std::vector<pt::idx_t> &depth = tm.depth;

  pt::idx_t ai_st_idx = ft_v.get_ai();
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
          pgt::id_or_t g = ft_v.get_end();
          res.push_back(pvst::Vertex::make_smothered(e, g, src, ft_v));
          //e.push_back(be.get_src());
          continue;
        }
      }

    }
  }

  // for (auto &v_idx : e) {
  //   std::cerr << fn_name << " found [1(b)]: " << st.get_vertex(sl_st_idx).g_v_id()
  //             << "~>" << st.get_vertex(v_idx).g_v_id() << "\n";
  // }
}

} // namespace g

namespace misc {
void trunk_misc_zi(const pst::Tree &st, const pvst::Vertex &ft_v, const ptu::tree_meta &tm) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  const std::vector<pt::idx_t> &depth = tm.depth;
  std::vector<pt::idx_t> e;

  pt::idx_t ai = ft_v.get_ai();
  pt::idx_t zi = ft_v.get_zi();

  pt::idx_t hi_trunk = zi;
  for (auto tgt : st.get_obe_tgt_v_idxs(zi)) {

    if (tgt == ai) {
      continue;
    }

    if (depth[tgt] < depth[hi_trunk]) {
      hi_trunk = tgt;
    }
  }

  if (hi_trunk == zi) {
    return;
  }

  for (auto src_v_idx : st.get_ibe_src_v_idxs(hi_trunk)) {

  if (depth[src_v_idx] < depth[zi]) {
    std::vector v {zi, src_v_idx};
    pt::idx_t lca = ptu::find_lca(tm, v);

    std::cerr << fn_name << " 2(a) (i): " << st.get_vertex(lca).g_v_id() << "\n";
    }

  }

}

void trunk_misc_ai(const pst::Tree &st, const pvst::Vertex &ft_v, const ptu::tree_meta &tm) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  const std::vector<pt::idx_t> &depth = tm.depth;
  std::vector<pt::idx_t> e;

  pt::idx_t ai = ft_v.get_ai();
  pt::idx_t zi = ft_v.get_zi();

  if(!(ai == 764 && zi == 765)) {
    return;
  }

  std::cerr << ft_v.as_str() << "\t" << " ai: " << ai << " zi: " << zi << "\n";

  for (auto src : st.get_ibe_src_v_idxs(ai)) {
    std::cerr << fn_name << " src: " << st.get_vertex(src).g_v_id() << "\n";
  }
}
} // namespace misc

namespace s {

void trunk(const pst::Tree &st, const pvst::Vertex &ft_v,
           const ptu::tree_meta &tm, std::vector<pvst::Vertex> &res) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

    const std::vector<pt::idx_t> &depth = tm.depth;
    //std::vector<pt::idx_t> e;

  pt::idx_t zi_st_idx = ft_v.get_zi();
  pt::idx_t sl_st_idx = ft_v.get_sl_st_idx();

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
      pgt::id_or_t w =
          (st.get_vertex(src).type() == pst::v_type_e::l)
              ? pgt::id_or_t{st.get_vertex(src).g_v_id(), pgt::or_e::reverse}
              : pgt::id_or_t{st.get_vertex(src).g_v_id(), pgt::or_e::forward};
      pgt::id_or_t s = ft_v.get_start();
      res.push_back(pvst::Vertex::make_smothered(s, w, src, ft_v));
      //e.push_back(src);
    }
  }


  // for (auto &v_idx : e) {
  //   std::cerr << fn_name << " found [2]: " << st.get_vertex(sl_st_idx).g_v_id()
  //             << "~>" << st.get_vertex(v_idx).g_v_id() << "\n";
  // }
}

void branch(const pst::Tree &st, const pvst::Vertex &ft_v,
            const ptu::tree_meta &tm, std::vector<pvst::Vertex> &res) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  //std::vector<pt::idx_t> e;

  pt::idx_t zi_st_idx = ft_v.get_zi();
  pt::idx_t sl_st_idx = ft_v.get_sl_st_idx();

  std::set<pt::idx_t> srcs = st.get_ibe_src_v_idxs(sl_st_idx);

  if (srcs.empty()) {
    return; // no srcs
  }

  for (pt::idx_t src_v_idx : srcs) {
    for (auto src_ : st.get_ibe_src_v_idxs(src_v_idx)){
      pgt::id_or_t w =
          (st.get_vertex(src_).type() == pst::v_type_e::l)
              ? pgt::id_or_t{st.get_vertex(src_).g_v_id(), pgt::or_e::reverse}
              : pgt::id_or_t{st.get_vertex(src_).g_v_id(), pgt::or_e::forward};
      pgt::id_or_t s = ft_v.get_start();
      res.push_back(pvst::Vertex::make_smothered(s, w, src_, ft_v));
      //e.push_back(src_);
    }

    for (auto tgt_ : st.get_obe_tgt_v_idxs(src_v_idx)){
      if (tgt_ == zi_st_idx) {
        continue; // skip the trunk
      }
      pgt::id_or_t w =
          (st.get_vertex(tgt_).type() == pst::v_type_e::l)
              ? pgt::id_or_t{st.get_vertex(tgt_).g_v_id(), pgt::or_e::reverse}
              : pgt::id_or_t{st.get_vertex(tgt_).g_v_id(), pgt::or_e::forward};
      pgt::id_or_t s = ft_v.get_start();
      res.push_back(pvst::Vertex::make_smothered(s, w, tgt_, ft_v));
      //e.push_back(tgt_);
    }
  }

  // for (auto &v_idx : e) {
  //   std::cerr << fn_name << " found [2(b)]: " << st.get_vertex(sl_st_idx).g_v_id() << "~>" << st.get_vertex(v_idx).g_v_id() << "\n";
  // }
}

} // namespace s

void add_smothered(const pst::Tree &st, pvtr::Tree<pvst::Vertex> &vst,
                   const ptu::tree_meta &tm, const std::vector<fl_sls> &al_smo) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  for (const fl_sls &smo : al_smo) {

    const pvst::Vertex &ft_v = vst.get_vertex(smo.fl_v_idx);

    if (ft_v.get_type() != pvst::vt_e::slubble) {
      std::cerr << fn_name << " found smothered for non-flubble vertex: " << smo.fl_v_idx << "\n";
      std::cerr << "flubble vertex: " << ft_v.as_str() << "\n";
      continue; // skip non-flubble vertices
    }

    for (pvst::Vertex g_adj_v : smo.g_adj) {
      std::cerr << fn_name << " adding smothered vertex: " << g_adj_v.as_str() << "\n";
      // use a ref and move?
      pt::idx_t g_adj_v_idx = vst.add_vertex(g_adj_v);
      vst.add_edge(smo.fl_v_idx, g_adj_v_idx);
      //vst.get_vertex(g_adj_v_idx).set_type(pvst::vt_e::smothered);
    }

    for (pvst::Vertex s_adj_v : smo.s_adj) {
      std::cerr << fn_name << " adding smothered vertex: " << s_adj_v.as_str()
                << "\n";
      // use a ref and move?
      pt::idx_t s_adj_v_idx = vst.add_vertex(s_adj_v);
      vst.add_edge(smo.fl_v_idx, s_adj_v_idx);
      // vst.get_vertex(g_adj_v_idx).set_type(pvst::vt_e::smothered);
    }
  }
}

void find_smothered(const pst::Tree &st, pvtr::Tree<pvst::Vertex> &ft, const ptu::tree_meta &tm) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  std::vector<fl_sls> all_smo;

  for (pt::idx_t ft_v_idx{}; ft_v_idx < ft.vtx_count(); ft_v_idx++) {
     const pvst::Vertex &ft_v = ft.get_vertex(ft_v_idx);

     if (ft_v.get_type() == pvst::vt_e::flubble) {
       misc::trunk_misc_ai(st, ft_v, tm);
       misc::trunk_misc_zi(st, ft_v, tm);
     }

     if (ft_v.get_type() != pvst::vt_e::slubble) {
       continue;
     }

     fl_sls smo {ft_v_idx};
     switch (ft_v.get_sl_type()) {
     case pvst::sl_type_e::ai_trunk:
       g::trunk(st, ft_v, tm, smo.g_adj);
       break;
     case pvst::sl_type_e::ai_branch:
       g::branch(st, ft_v, tm, smo.g_adj);
       break;
     case pvst::sl_type_e::zi_trunk:
       s::trunk(st, ft_v, tm, smo.s_adj);
       break;
     case pvst::sl_type_e::zi_branch:
       s::branch(st, ft_v, tm, smo.s_adj);
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
