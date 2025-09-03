#include "./tiny.hpp"

namespace povu::tiny {

/**
 * Y = { C(zi) \ zx }
 */
std::vector<pt::idx_t> find_Y(const pst::Tree &st, pt::idx_t zi) {
  std::vector<pt::idx_t> Y;
  for (auto c_e : st.get_child_edge_idxs(zi)) {
    const pst::Edge &c_e_ref = st.get_tree_edge(c_e);
    if (c_e_ref.get_color() == pgt::color_e::black) {
      continue;
    }

    Y.push_back(c_e_ref.get_child_v_idx());
  }

  return Y;
}

bool branches(const pst::Tree &st, const ptu::tree_meta &tm, pt::idx_t ai, pt::idx_t zi) {
  const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  // all brackets are to ai
  auto has_be_to_ai = [&](pt::idx_t c) -> bool {
    for (pt::idx_t be_idx : tm.get_brackets(c)) {
      if (ai == st.get_be(be_idx).get_tgt()) {
        return true;
      }
    }

    for (pt::idx_t tgt_v_idx : st.get_obe_idxs(c)) {
      if (ai == tgt_v_idx) {
        return true;
      }
    }

    return false;
  };

  // the vertex c has exactly one descendant
  auto has_one_descendants = [&](pt::idx_t v_idx) -> bool {
    pt::idx_t x = st.get_vertex(v_idx).post_order() - st.get_vertex(v_idx).pre_order();
    return x == 3;
  };

  std::vector<pt::idx_t> Y = find_Y(st, zi);
  for (pt::idx_t c_v_idx : Y) {
    if (!has_one_descendants(c_v_idx) || !has_be_to_ai(c_v_idx)) {
      return false;
    }
  }

  return true;
}

bool trunk(const pst::Tree &st, pt::idx_t ai, pt::idx_t zi) {
  const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  for (auto c_e_idx : st.get_child_edge_idxs(zi)) {
    const pst::Edge &e = st.get_tree_edge(c_e_idx);
    if (e.get_color() == pgt::color_e::black) {
      continue;
    }

    return false;
  }

  // TODO: [A] investigate
  if (st.get_obe_tgt_v_idxs(zi).empty()) {
    return false;
  }

#ifdef DEBUG
  assert(st.get_obe_tgt_v_idxs(zi).size() == 1);
  assert(*(st.get_obe_tgt_v_idxs(zi).begin()) == ai);
#endif

  return false;
}

void find_tiny(const pst::Tree &st, pvst::Tree &ft, const ptu::tree_meta &tm) {
  const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  for (pt::idx_t ft_v_idx{}; ft_v_idx < ft.vtx_count(); ft_v_idx++) {

    if (!ft.is_leaf(ft_v_idx)) {
      continue;
    }

    pvst::VertexBase &pvst_v = ft.get_vertex_mut(ft_v_idx);

    if (pvst_v.get_fam() != pvst::vt_e::flubble) {
      continue;
    }

    pvst::Flubble &ft_v = static_cast<pvst::Flubble&>(pvst_v);

    pt::idx_t ai = ft_v.get_ai();
    pt::idx_t zi = ft_v.get_zi();

     if (!(zi - ai == 1 || zi - ai == 3)) {
      continue;
    }

    if (trunk(st, ai, zi) || branches(st, tm, ai, zi)) {
      ft_v.set_type(pvst::vt_e::tiny);
    }
  }
}

} // namespace povu::tiny
