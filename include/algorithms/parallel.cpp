#include "./parallel.hpp"

namespace povu::parallel {

bool inspect_trunk(const pst::Tree &st, pt::idx_t ai, pt::idx_t zi) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  //bool dbg = (ai == 1338 && zi == 1343) ? true : false;

  // go up from zi to ai and ensure no branching vertices
  auto cond_a = [&]() -> bool {
    pt::idx_t v_idx = zi;
    while (v_idx != ai) {
      if (st.get_child_count(v_idx) > 1) {
        return false;
      }
      v_idx = st.get_parent_v_idx(v_idx);
    }

    return true;
  };

  auto cond_b = [&]() -> bool {
    pt::idx_t branching_vtx {pc::INVALID_IDX};
    pt::idx_t v_idx = zi;

    while (v_idx != ai) {
      if (st.get_child_count(v_idx) > 1) {
        if (branching_vtx == pc::INVALID_IDX){
          branching_vtx = v_idx;
        } else {
          // more than one branching vertex
          return false;
        }

      }
      v_idx = st.get_parent_v_idx(v_idx);
    }

    if (branching_vtx == pc::INVALID_IDX) { // no branching vertex found
      // if (dbg) {
      //   std::cerr << "No branching vertex found from " << st.get_vertex(zi).g_v_id() << " to " << st.get_vertex(ai).g_v_id() << "\n";
      // }
      return false;
    }
    //else {
      // if (dbg) {
      //   std::cerr << "Branching vertex found: " << st.get_vertex(branching_vtx).g_v_id() << "\n";
      // }
    //}

    for (pt::idx_t c_v_idx : st.get_children(branching_vtx)) {
      // if(dbg ) {
      //   std::cerr << "c " << c_v_idx << " id "
      //             << st.get_vertex(c_v_idx).g_v_id() << " post order"
      //             << st.get_vertex(c_v_idx).post_order() << " pre order "<<
      //                    st.get_vertex(c_v_idx).pre_order()
      //             << "\n";
      // }
      if (c_v_idx > zi && (st.get_vertex(c_v_idx).post_order() -  st.get_vertex(c_v_idx).pre_order() == 3)) {
        return true;
      }
    }

    return false;
  };

  // if (dbg) {
  //   std::cerr << "[DEBUG] Inspecting trunk from " << st.get_vertex(zi).g_v_id() << " to " << st.get_vertex(ai).g_v_id() << "\n";
  //   std::cerr << "b " << cond_b() << "\n";
  //   exit(1);
  // }

   return cond_a() ||  cond_b();
 }


bool in_trunk(const pst::Tree &st, const pvst::Flubble &ft_v) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  pt::idx_t ai = ft_v.get_ai();
  pt::idx_t zi = ft_v.get_zi();

  // if (zi < 3) {
  //   return false;
  // }

  // condition i
  if (zi - ai <= 3 && st.get_ibe_idxs(ai).size() <= 1) {
    return false;
  }

  // condition iii
  if (st.get_child_count(zi) != 1) {
    return false;
  }

  // condition iv
  if (!inspect_trunk(st, ai, zi)) {
    return false;
  }

  // with ai
  {
    pt::idx_t be_count {};
    // count OBE(zi) \ capping be
    // for (auto be_idx : st.get_obe_idxs(zi)) {
    //   const pst::BackEdge &be = st.get_be(be_idx);
    //    if (be.type() != pst::be_type_e::back_edge) {
    //     continue;
    //   }

    //   be_count++;
    // }

    // if (be_count < 2) {
    //   return false;
    // }

    //be_count = 0; // reset be_count
    // count IBE(ai) \ capping be
    for (auto be_idx : st.get_ibe_idxs(ai)) {
      const pst::BackEdge &be = st.get_be(be_idx);
      if (be.type() != pst::be_type_e::back_edge) {
        continue;
      }

      be_count++;
    }


    if (2 * be_count >= ((zi - ai) - 3)) {
      return true;
    }
  }

  // with zi
  {
    pt::idx_t be_count{};
    // count IBE(ai) \ capping be
    for (auto be_idx : st.get_ibe_idxs(ai)) {
      const pst::BackEdge &be = st.get_be(be_idx);
      if (be.type() != pst::be_type_e::back_edge) {
        continue;
      }

      be_count++;
    }

    if (be_count != 0) {
      return false;
    }

    be_count = 0; // reset be_count
    // count OBE(zi) \ capping be
    for (auto be_idx : st.get_obe_idxs(zi)) {
      const pst::BackEdge &be = st.get_be(be_idx);
      if (be.type() != pst::be_type_e::back_edge) {
        continue;
      }

      be_count++;
    }

    if (2 * be_count >= (zi - ai) - 3) {
      return true;
    }
  }

  return false;
}


bool in_branch(const pst::Tree &st, const ptu::tree_meta &tm, const pvst::Flubble &ft_v) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  pt::idx_t ai = ft_v.get_ai();
  pt::idx_t zi = ft_v.get_zi();

  // if (zi < 3) {
  //   return false;
  // }

  if (zi - ai != 1) {
    return false;
  }

  // ensure only 1 gray child from zi
  pt::idx_t c_idx{pc::INVALID_IDX};
  for (pt::idx_t e_idx : st.get_child_edge_idxs(zi)) {
    const pst::Edge &e = st.get_tree_edge(e_idx);
    if (e.get_color() == pgt::color_e::black) {
      continue;
    }

    // Y != 1
    if (c_idx != pc::INVALID_IDX) {
      return false;
    }

    c_idx = e.get_child_v_idx();
  }

  pt::idx_t br_count = tm.get_brackets(c_idx).size();
  pt::idx_t ch_obe_count = st.get_obe_idxs(c_idx).size();

  if (br_count <= 2) {
    return false;
  }

  // with ai
  {
    pt::idx_t be_count{};
    // count IBE(ai) \ capping be
    for (auto be_idx : st.get_ibe_idxs(ai)) {
      const pst::BackEdge &be = st.get_be(be_idx);
      if (be.type() != pst::be_type_e::back_edge) {
        continue;
      }

      be_count++;
    }

    if (be_count >= br_count + ch_obe_count) {
      return true;
    }
  }

  // with zi
  {
    pt::idx_t be_count{};
    be_count = st.get_obe_idxs(zi).size(); // reset be_count

    if (be_count >= br_count + ch_obe_count) {
      return true;
    }
  }

  return false;
}


void find_parallel(const pst::Tree &st, pvtr::Tree &ft, const ptu::tree_meta &tm) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  for (pt::idx_t ft_v_idx{}; ft_v_idx < ft.vtx_count(); ft_v_idx++) {
    if (!ft.is_leaf(ft_v_idx)) {
      continue;
    }

    pvst::VertexBase &pvst_v = ft.get_vertex_mut(ft_v_idx);

    if (pvst_v.get_type() != pvst::vt_e::flubble) {
      continue;
    }

    pvst::Flubble &ft_v = static_cast<pvst::Flubble &>(pvst_v);

    if (in_branch(st, tm, ft_v) || in_trunk(st, ft_v)) {
      ft_v.set_type(pvst::vt_e::parallel);
      // std::cerr << fn_name << ": Found parallel flubble at " << ft_v.as_str() << "\n";
    }
  }
}

} // namespace povu::parallel
