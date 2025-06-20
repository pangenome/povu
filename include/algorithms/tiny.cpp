#include "./tiny.hpp"

namespace povu::tiny {

bool handle_branches(const pst::Tree &st, const ptu::tree_meta &tm,
                     const pvst::Flubble &ft_v, std::vector<pt::idx_t> Y) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  pt::idx_t ai = ft_v.get_ai();

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

  for (pt::idx_t c : Y) {
    pt::idx_t x = st.get_vertex(c).post_order() - st.get_vertex(c).pre_order();

    // std::cerr << fn_name << " " << st.get_vertex(c).g_v_id() << " x: " << x << "\n";

    if (x != 3)  {
      return false;
    }

    if (!has_be_to_ai(c)) {
      return false;
    }
  }

  return true;
}


void find_tiny(const pst::Tree &st, pvtr::Tree &ft, const ptu::tree_meta &tm) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  for (pt::idx_t ft_v_idx{}; ft_v_idx < ft.vtx_count(); ft_v_idx++) {
    pvst::VertexBase &pvst_v = ft.get_vertex_mut(ft_v_idx);

    if (pvst_v.get_type() != pvst::vt_e::flubble) {
      continue;
    }

    pvst::Flubble &ft_v = static_cast<pvst::Flubble&>(pvst_v);

    pt::idx_t ai = ft_v.get_ai();
    pt::idx_t zi = ft_v.get_zi();

    if (!(zi - ai == 1 || zi-ai == 3)) {
      continue;
    }

    std::vector<pt::idx_t> Y;
    for (auto c_e : st.get_child_edge_idxs(zi)){
      const pst::Edge &c_e_ref = st.get_tree_edge(c_e);
      if (c_e_ref.get_color() == pgt::color_e::black) {
        continue;
      }

      Y.push_back(c_e_ref.get_child_v_idx());
    }

    if (Y.empty()) {
#ifdef DEBUG
      assert(st.get_obe_tgt_v_idxs(zi).size() == 1);
      assert(*(st.get_obe_tgt_v_idxs(zi).begin()) == ai);
#endif
      ft_v.set_type(pvst::vt_e::tiny);
      //std::cerr << fn_name << ": Found tiny flubble at " << ft_v.as_str() << "\n";
      continue;
    }

    if (handle_branches(st, tm, ft_v, Y)) {
      ft_v.set_type(pvst::vt_e::tiny);
      //std::cerr << fn_name << ": Found tiny flubble at " << ft_v.as_str() << "\n";
    }
  }
}

  } // namespace povu::parallel
 
