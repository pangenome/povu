#include "./misc.hpp"

namespace povu::misc {
struct fl_sls {
  pt::idx_t cn_v_idx;
  std::vector<pvst::Smothered> g_adj;
  std::vector<pvst::Smothered> s_adj;

  pt::idx_t size() const { return g_adj.size() + s_adj.size(); }

  // Constructor for fl_sls
  fl_sls(pt::idx_t cn_v_idx_)
      : cn_v_idx(cn_v_idx_), g_adj(std::vector<pvst::Smothered>{}),
        s_adj(std::vector<pvst::Smothered>{}) {}
};

void trunk_misc_zi(const pst::Tree &st, const pvst::Flubble &ft_v, const ptu::tree_meta &tm) {
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

void trunk_misc_ai(const pst::Tree &st, const pvst::Flubble &ft_v, const ptu::tree_meta &tm) {
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

void find_misc(const pst::Tree &st, pvtr::Tree &pvst, const ptu::tree_meta &tm) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  std::vector<fl_sls> all_misc;

  for (pt::idx_t ft_v_idx{}; ft_v_idx < pvst.vtx_count(); ft_v_idx++) {
     const pvst::VertexBase &pvst_v = pvst.get_vertex(ft_v_idx);

     if (pvst_v.get_type() != pvst::vt_e::flubble) {
       continue;
     }

     const pvst::Flubble &fl_v = static_cast<const pvst::Flubble &>(pvst_v);

     misc::trunk_misc_ai(st, fl_v, tm);
     misc::trunk_misc_zi(st, fl_v, tm);
  }
}
} // namespace misc
