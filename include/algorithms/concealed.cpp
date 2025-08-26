#include "./concealed.hpp"


namespace povu::concealed {

/*
  -----------------
  Types
  -----------------
*/
// slubbles in flubble or flubble slubbles
struct fl_sls {
  pt::idx_t fl_v_idx;
  std::vector<pvst::Concealed> ii_adj;
  std::vector<pvst::Concealed> ji_adj;
  pt::idx_t m;
  pt::idx_t n;

  pt::idx_t size() const {
    return ii_adj.size() + ji_adj.size();
  }
};


struct mn_t {
  pt::idx_t m;
  pt::idx_t n;
};


struct src_lca_t {
  pt::idx_t be_src; // be src vertex index
  pt::idx_t lca; // ℓ
};


const src_lca_t INVALID_SRC_LCA{pc::INVALID_IDX, pc::INVALID_IDX};

/*
  -----------------
  Utils
  -----------------
*/
/**
 * @brief check if a is a descendant of d
 * @param st the pst tree
 * @param a the potential ancestor vertex index
 * @param d the potential descendant vertex index
 * @return true if a is an ancestor of d, false otherwise
 */
bool is_desc(const pst::Tree &st, pt::idx_t a, pt::idx_t d) {
  return st.get_vertex(a).pre_order() < st.get_vertex(d).pre_order() &&
         st.get_vertex(a).post_order() > st.get_vertex(d).post_order();
}

pvst::Concealed gen_ai_slubble(const pst::Tree &st, pt::idx_t ai_st_v_idx,
                               src_lca_t tb, pvst::cl_e t, pt::idx_t fl_v_idx) {
  const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  auto [be_src_v_idx, sl_st_idx]= tb;

  // find a
  pst::Vertex st_v = st.get_vertex(ai_st_v_idx);
  pt::id_t fl_id = st_v.g_v_id();
  pgt::or_e fl_o = (st_v.type() == pgt::v_type_e::r) ? pgt::or_e::forward : pgt::or_e::reverse;
  pgt::id_or_t fl_boundary = pgt::id_or_t{fl_id, fl_o};

  // compute g
  pt::idx_t sl_id = st.get_vertex(sl_st_idx).g_v_id();
  pgt::or_e sl_o;
  if (t == pvst::cl_e::ai_trunk) {
    if (st.get_parent_edge(sl_st_idx).get_color() == pgt::color_e::black) {
      sl_o = st.get_vertex(sl_st_idx).type() == pgt::v_type_e::r ? pgt::or_e::forward : pgt::or_e::reverse;
    }
    else { // has a black child
      sl_o = st.get_vertex(sl_st_idx).type() == pgt::v_type_e::r ? pgt::or_e::reverse : pgt::or_e::forward;
    }
  }
  else if (t == pvst::cl_e::ai_branch) {
    if (st.get_parent_edge(sl_st_idx).get_color() == pgt::color_e::black) {
      sl_o = st.get_vertex(sl_st_idx).type() == pgt::v_type_e::r ? pgt::or_e::reverse : pgt::or_e::forward;
    } else { // has a black child
      sl_o = st.get_vertex(sl_st_idx).type() == pgt::v_type_e::r ? pgt::or_e::forward : pgt::or_e::reverse;
    }
  }
  else {
    std::string err_msg = pv_cmp::format("{} called with invalid slubble type", fn_name);
    perror(err_msg.c_str());
    exit(1);
  }

  //
  pgt::id_or_t sl_boundary = pgt::id_or_t{sl_id, sl_o};
  pgt::id_or_t a;
  pgt::id_or_t g;
  if (sl_o == pgt::or_e::reverse && fl_o == pgt::or_e::reverse) {
    a = pgt::id_or_t{sl_id, pgt::or_e::forward};
    g = pgt::id_or_t{fl_id, pgt::or_e::forward};
  }
  else {
    a = fl_boundary;
    g = sl_boundary;
  }

  pvst::bounds_t bounds;
  if (t == pvst::cl_e::ai_trunk) {
    bounds = is_desc(st, ai_st_v_idx, be_src_v_idx)
                 ? pvst::bounds_t{ai_st_v_idx, be_src_v_idx}
                 : pvst::bounds_t{be_src_v_idx, ai_st_v_idx};
  }
  else {
    bounds = { sl_st_idx, pc::INVALID_IDX };
  }

  return pvst::Concealed::create(a, g, bounds, fl_v_idx, t, sl_st_idx, pvst::rt_e::e2s);
}

pvst::Concealed gen_zi_slubble(const pst::Tree &st, pt::idx_t zi_st_v_idx,
                            pt::idx_t sl_st_idx, pvst::cl_e t,
                            pt::idx_t fl_v_idx) {
  const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  // find z
  pst::Vertex zi_st_v = st.get_vertex(zi_st_v_idx);
  pt::id_t fl_id = zi_st_v.g_v_id();
  pgt::or_e fl_o = (zi_st_v.type() == pgt::v_type_e::l) ? pgt::or_e::forward : pgt::or_e::reverse;
  pgt::id_or_t fl_boundary = pgt::id_or_t{fl_id, fl_o};

  // compute s
  pt::idx_t sl_id = st.get_vertex(sl_st_idx).g_v_id();
  pgt::or_e sl_o;
  if (st.get_parent_edge(sl_st_idx).get_color() == pgt::color_e::black) {
    sl_o = st.get_vertex(sl_st_idx).type() == pgt::v_type_e::r ?  pgt::or_e::forward : pgt::or_e::reverse;
  }
  else { // has a black child
    sl_o = st.get_vertex(sl_st_idx).type() == pgt::v_type_e::r ? pgt::or_e::reverse : pgt::or_e::forward;
  }

  //
  pgt::id_or_t sl_boundary = pgt::id_or_t{sl_id, sl_o};
  pgt::id_or_t s;
  pgt::id_or_t z;
  if (sl_o == pgt::or_e::reverse && fl_o == pgt::or_e::reverse) {
    s = pgt::id_or_t{sl_id, pgt::or_e::forward};
    z = pgt::id_or_t{fl_id, pgt::or_e::forward};
  }
  else {
    s = sl_boundary;
    z = fl_boundary;
  }

  pvst::bounds_t bounds = is_desc(st, zi_st_v_idx, sl_st_idx)
    ? pvst::bounds_t{zi_st_v_idx, sl_st_idx}
    : pvst::bounds_t{sl_st_idx, zi_st_v_idx};

  return pvst::Concealed::create(z, s, bounds, fl_v_idx, t, sl_st_idx, pvst::rt_e::s2e);
}


// not inclusive
bool is_btwn(const pst::Tree &st, pt::idx_t v_idx, pt::idx_t upper, pt::idx_t lower) {
  return is_desc(st, upper, v_idx) && is_desc(st, v_idx, lower);
}


pt::idx_t compute_m(const pst::Tree &st, const ptu::tree_meta &tm,
                    pt::idx_t ii_v_idx, pt::idx_t ji_v_idx) {
   const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};
  if (st.get_ibe_src_v_idxs(ii_v_idx).size() == 0) {
    return ii_v_idx;
  }

  const std::vector<pt::idx_t> &height = tm.depth;

  struct v_idx_height_t {
    pt::idx_t v_idx;
    pt::idx_t height;
  };

  v_idx_height_t max_height{ pc::INVALID_IDX, 0 };

  for (pt::idx_t be_idx : st.get_ibe_idxs(ii_v_idx)) {
    if (st.get_be(be_idx).type() != pst::be_type_e::back_edge) {
      continue;
    }

    pt::idx_t src_v_idx = st.get_be(be_idx).get_src();
    std::vector<pt::idx_t> p {src_v_idx, ji_v_idx};
    pt::idx_t l = find_lca(tm, p);

    if (height[l] >= height[ji_v_idx]) {
      // std::cerr << fn_name << "hight[src] " << height[src_v_idx] << " h[ji] " << height[ji_v_idx] << "\n";
      // std::cerr << fn_name << " skipping " << src_v_idx << "\n";
      continue;
    }

    if (height[l] > max_height.height) {
      // TODO: [c] remove type cast
      max_height = {(pt::idx_t)l, height[l]};
    }
  }

  return  max_height.v_idx == pc::INVALID_IDX ? ii_v_idx : max_height.v_idx;
}


pt::idx_t compute_n(const pst::Tree &st, const ptu::tree_meta &tm,
                    pt::idx_t ii_v_idx, pt::idx_t ji_v_idx) {

  pt::idx_t count{0};
  // count normal OBE
  for (pt::idx_t be_idx : st.get_obe_idxs(ji_v_idx)) {
    if (st.get_be(be_idx).type() == pst::be_type_e::back_edge) {
      count++;
    }
  }

  if (count == 0) {
    return ji_v_idx;
  }

  //std::cerr << "OBE cnt " << st.get_obe_tgt_v_idxs(ji_v_idx).size() << "\n";

  const std::vector<pt::idx_t> &depth = tm.depth;


  pt::idx_t lowest { ii_v_idx };


  for (pt::idx_t be_idx : st.get_obe_idxs(ji_v_idx)) {
    const pst::BackEdge &be = st.get_be(be_idx);
    if (be.type() != pst::be_type_e::back_edge) {
      continue;
    }

    pt::idx_t tgt_v_idx = be.get_tgt();

    //if (ii_v_idx == tgt_v_idx) {
    //  continue;
    //}

    if (depth[tgt_v_idx] > depth[lowest]) {
      lowest = tgt_v_idx;
    }
  }

  if (lowest == ii_v_idx) {
    return ji_v_idx;
  }

  return lowest;
}


mn_t get_mn(const pst::Tree &st, const ptu::tree_meta &tm, pt::idx_t ii_v_idx,
            pt::idx_t ji_v_idx) {
  // when m is invalid (i) of ii trunk be is false

  pt::idx_t m = compute_m(st, tm, ii_v_idx, ji_v_idx);
  pt::idx_t n = compute_n(st, tm, ii_v_idx, ji_v_idx);

  return {m, n};
}


/**
 * check if a flubble can not contain a slubble
 * @param tm the tree meta data
 * @param ii_v_idx
 * @param ji_v_idx
 * @return true if the flubble can not contain the slubble, false otherwise
 */
bool can_contain(const pst::Tree &st,  const ptu::tree_meta &tm,
                 pt::idx_t ii_v_idx, pt::idx_t ji_v_idx, pt::idx_t m,
                 pt::idx_t n) {
  const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  const std::vector<pt::idx_t> &depth = tm.depth;
  const std::vector<pt::idx_t> &lo = tm.lo;

  auto count_ibe_ii = [&]() {
    pt::idx_t count {};
    for (pt::idx_t be_idx : st.get_ibe_idxs(ii_v_idx)) {
      const pst::BackEdge &be = st.get_be(be_idx);
      if (be.type() == pst::be_type_e::back_edge && be.get_src() != ji_v_idx) {
        count++;
      }
    }
    return count;
  };

  auto count_ibe_ji = [&]() {
    pt::idx_t count{};
    for (pt::idx_t be_idx : st.get_ibe_idxs(ji_v_idx)) {
      if (st.get_be(be_idx).type() == pst::be_type_e::back_edge) {
        count++;
      }
    }
    return count;
  };

  // (i)
  if (depth[lo[ji_v_idx]] < depth[ii_v_idx]) {
    return false;
  }

  // (ii)
  if (m == ii_v_idx && n == ji_v_idx && count_ibe_ii() < 2 &&
      count_ibe_ji() == 0 && st.get_children(ji_v_idx).size() < 3) {
    return false;
  }

  return true;
}


namespace ai {
src_lca_t ai_trunk(const pst::Tree &st, const ptu::tree_meta &tm, pt::idx_t m,
                   pt::idx_t n, pt::idx_t ai, pt::idx_t zi) {
  const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  const std::vector<pt::idx_t> &depth = tm.depth;

  //bool dbg = (ai == 1664 && zi == 1669) ? true : false;

  if (n == pc::INVALID_IDX || m == pc::INVALID_IDX || depth[m] > depth[n]) {
    return INVALID_SRC_LCA;
  }

  // cond ii brackets of ell
  auto ell_brackets = [&](pt::idx_t l) -> bool {
    for (pt::idx_t be_idx : tm.get_brackets(l)) {
      const pst::BackEdge &be = st.get_be(be_idx);
      //pt::idx_t br_src_v_idx = be.get_src();
      pt::idx_t br_tgt_v_idx = be.get_tgt();

      bool is_tgt_ai_or_above = depth[br_tgt_v_idx] <= depth[ai];
      //bool is_src_branch_vtx = is_desc(st, zi, br_src_v_idx);
      //bool is_src_btwn_m_zi = is_btwn(st, br_src_v_idx, m, zi);

      // if(dbg) {
      //   std::cerr << fn_name << " tgt <= ai " << is_tgt_ai_or_above <<  " tgt " << br_tgt_v_idx << " src" << br_src_v_idx << "\n";
      // }

      if (!is_tgt_ai_or_above ) {
        return false;
      }
      // if (!is_tgt_ai_or_above && (is_src_branch_vtx ||is_src_btwn_m_zi)){
      //   return false;
      // }
   }
    return true;
  };

   auto cond_iii = [&](pt::idx_t ell) -> bool {
    return !st.get_obe_idxs(ell).empty() || st.get_child_count(ell) > 1;
  };

  // cond ii and iii
  auto not_cond_ii_iii = [&](src_lca_t x) -> bool {
    return !ell_brackets(x.lca) && cond_iii(x.lca);
  };

  //bool dbg = (ai == 1772) ? true : false;

  // pt::idx_t v = zi;
  // std::cerr << fn_name << " trunk\n";
  // while(dbg && v >= ai) {
  //   std::cerr << "\t " << st.get_vertex(v).g_v_id() << ", ";
  //   // print children
  //   std::cerr << "[ ";
  //   for (pt::idx_t c_v_idx : st.get_children(v)) {
  //     std::cerr << st.get_vertex(c_v_idx).g_v_id() << ", ";
  //   }
  //   std::cerr << "]\n";
  //   v = st.get_parent_v_idx(v);
  // }

  // cond i
  // ------
  std::vector<src_lca_t> src_lca_vec;
  for (pt::idx_t be_idx : st.get_ibe_idxs(ai)) {
    const pst::BackEdge &be = st.get_be(be_idx);

    if (be.type() != pst::be_type_e::back_edge) {
      continue;
    }

    pt::idx_t be_src_v_idx = be.get_src();
    std::vector<pt::idx_t> p{be_src_v_idx, zi};
    pt::idx_t l = find_lca(tm, p);



    // if (dbg) {
    //   std::cerr << " l " << l << " ai " << ai << "\n";

    //   std::cerr << fn_name << " ell br " << ell_brackets(l) << "\n";

    //   pt::idx_t v = l;
    //   //std::cerr << fn_name << " trunk\n";

    //   // all brackets are from below zi

    //   while(dbg && v >= ai) {
    //     std::cerr << "\t " << st.get_vertex(v).g_v_id() << ", ";
    //     // print children
    //     std::cerr << "[ ";
    //     for (pt::idx_t c_v_idx : st.get_children(v)) {
    //       std::cerr << st.get_vertex(c_v_idx).g_v_id() << ", ";
    //     }
    //     std::cerr << "] \t";

    //     // print IBE
    //     std::cerr << "IBE: (";
    //     for(auto src : st.get_ibe_src_v_idxs(v)){
    //       std::cerr << st.get_vertex(src).g_v_id() << ", ";
    //     }
    //     std::cerr << ")\n";
    //     v = st.get_parent_v_idx(v);
    //   }
    // }




    if (depth[l] <= depth[m]) {
      // if (dbg) {
      //   std::cerr << fn_name << " " << l << " " << st.get_vertex(l).g_v_id() << "\n";
      // }

      src_lca_vec.push_back({be_src_v_idx, l});
    }
  }

  // remove elements that do not meet condition ii and iii
  pv_cmp::erase_if(src_lca_vec, not_cond_ii_iii);

  // erase_if is C++20
  // remove elements that do not meet condition iv and v
  // std::erase_if(src_lca_vec, not_cond_ii_iii);

  // remove_if is C++11
  //   std::remove_if(...) moves the elements that don't match the predicate to the end.
  //   .erase(...) trims those elements off the container.
  //src_lca_vec.erase(std::remove_if(src_lca_vec.begin(), src_lca_vec.end(), not_cond_ii_iii), src_lca_vec.end());

  // sort by LCA depth (same as dfs num in this case) in ascending order
  std::sort(src_lca_vec.begin(), src_lca_vec.end(),
                        [&](src_lca_t a, src_lca_t b) { return a.lca < b.lca; });

  if (!src_lca_vec.empty()) {
    return src_lca_vec.back(); // return the last element (with max depth)
  }

  return INVALID_SRC_LCA; // no valid trunk concealed bubble found
}

/** get Ii brach backedges */
void ai_branches(const pst::Tree &st, const ptu::tree_meta &tm,
                 pt::idx_t ai, pt::idx_t zi,
                 std::vector<pt::idx_t> &bb) {

  // condition (i)
  if (st.get_child_count(zi) < 2) {
    return;
  }

  const std::vector<pt::idx_t> &height = tm.D;
  const std::vector<pt::idx_t> &lo = tm.lo;

  // children who meet condition (ii)
  std::vector<pt::idx_t> branches;
  for (pt::idx_t c_v_idx : st.get_children(zi)) {
    if (st.get_vertex(c_v_idx).hi() == lo[c_v_idx] && st.get_vertex(c_v_idx).hi() == ai) {
      std::vector<pt::idx_t> c_br = tm.get_brackets(c_v_idx);
      if (c_br.size() > 1) {
        branches.push_back(c_v_idx);
      }
    }
  }

  for (auto c_v_idx : branches) {
    // std::cerr << "child" << c_v_idx << "\n";
    std::vector<pt::idx_t> brackets = tm.get_brackets(c_v_idx);

    // case (ii) (a)
    if (brackets.size() < 2) {
      continue;
    }

    std::vector<pt::idx_t> br_srcs;
    for (auto be_idx : brackets) {
      const pst::BackEdge &be = st.get_be(be_idx);
      pt::idx_t src_v_idx = be.get_src();
      br_srcs.push_back(src_v_idx);
    }

    pt::idx_t d = ptu::find_lca(tm, br_srcs);
    if (tm.get_brackets(d).size() > 0) {
      bb.push_back(d);
    }
    else {
      // get the br srcs with min depth
      // v_idx, height
      // pt::idx_t min;
      std::pair<pt::idx_t, pt::idx_t> min_depth {pc::INVALID_IDX, pc::MAX_IDX};
      for (pt::idx_t be_idx : tm.get_brackets(d)) {
        const pst::BackEdge &be = st.get_be(be_idx);
        pt::idx_t br_src = be.get_src();
        if (height[br_src] < min_depth.second) {
          min_depth = {br_src, height[br_src]};
        }
      }

      if (min_depth.first != pc::INVALID_IDX) {
        bb.push_back(min_depth.first);
      }
    }
  }
}

void with_ai(const pst::Tree &st, const ptu::tree_meta &tm,
             std::vector<pvst::Concealed> &res, pt::idx_t m, pt::idx_t n,
             pt::idx_t ii_v_idx, pt::idx_t ji_v_idx, pt::idx_t fl_v_idx) {
  const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  src_lca_t tb = ai_trunk(st, tm, m, n, ii_v_idx, ji_v_idx);
  if (tb.lca != pc::INVALID_IDX) {
    pvst::Concealed sl = gen_ai_slubble(st, ii_v_idx, tb, pvst::cl_e::ai_trunk, fl_v_idx);
    res.push_back(sl);
    //std::cerr << "trunk sl: " << sl.as_str() << "\n";
  }

  std::vector<pt::idx_t> bb; // branch boundaries
  ai_branches(st, tm, ii_v_idx, ji_v_idx, bb);
  //std::cerr << "branch sls:\n";
  for (auto b : bb) {
    src_lca_t tb_ {pc::INVALID_IDX, b}; // TODO: [A] fix the inv idx
    pvst::Concealed sl = gen_ai_slubble(st, ii_v_idx, tb_, pvst::cl_e::ai_branch, fl_v_idx);
    res.push_back(sl);
    //std::cerr << sl.as_str() << ", ";
  }
  // std::cerr << "\n";
}
}; // namespace ai


namespace zi {
pt::idx_t override_ji_trunk(const pst::Tree &st, const ptu::tree_meta &tm,
                            pt::idx_t m, pt::idx_t n,
                            pt::idx_t ii_v_idx, pt::idx_t ji_v_idx) {

  const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  //std::cerr << pv_cmp::format("{}\n", fn_name);

  //auto [m, n] = mn;
  const std::vector<pt::idx_t> &height = tm.depth; // rename to depth

  if (height[m] > height[n]) {
    return pc::INVALID_IDX; // invalid
  }


  // std::cerr << fn_name << " m: " << m << " height[m] " << height[m] << " n " << n
  //           << " height[n]: " << height[n] << "\n";



    std::vector<pt::idx_t> case_ii_a_y;
    for (pt::idx_t c_v_idx : st.get_children(ji_v_idx)) {

      //std::cerr << fn_name << " c : " << c_v_idx << "\n";

      if (height[st.get_vertex(c_v_idx).hi()] < height[ii_v_idx]) {
        // when true the child is j_x
        //std::cerr << fn_name << " skip child " << c_v_idx << "\n";
        continue;
      }

      //for (pt::idx_t be_idx : tm.get_brackets(c_v_idx)) {
        //const pst::BackEdge &be = st.get_be(be_idx);
        //pt::idx_t src_v_idx = be.get_src();
        //pt::idx_t tgt_v_idx = be.get_tgt();

        //std::cerr << fn_name << " src: " << src_v_idx << " tgt: " << tgt_v_idx
        //         << "\n";
        //}

      if (tm.get_brackets(c_v_idx).size() == 1) { // case (ii)
        std::vector<pt::idx_t> br = tm.get_brackets(c_v_idx);
        pt::idx_t be_idx = br.front();
        const pst::BackEdge &be = st.get_be(be_idx);
        //pt::idx_t src_v_idx = be.get_src();
        pt::idx_t tgt_v_idx = be.get_tgt();

        if (height[m] < height[tgt_v_idx] && height[n] > height[tgt_v_idx]) {
          case_ii_a_y.push_back(tgt_v_idx); // case (ii) (a)
        }
      }
    }

    pt::idx_t min_v{n};
    // no branching path from source of only bracket to above j_i
    for (auto v_idx : case_ii_a_y) {
      bool valid{true};
      for (pt::idx_t be_idx : tm.get_brackets(v_idx)) {
        const pst::BackEdge &be = st.get_be(be_idx);
        pt::idx_t src_v_idx = be.get_src();
        pt::idx_t tgt_v_idx = be.get_tgt();

        if (height[src_v_idx] < height[ji_v_idx] ||
            height[tgt_v_idx] > height[ii_v_idx]) {
          valid = false;
          break;
        }
      }

      if (valid && height[v_idx] < height[min_v]) {
        min_v = v_idx;
      }
    }

    if (min_v == n) {
      return pc::INVALID_IDX;
    } else {
      return min_v;
    }
}

pt::idx_t ji_trunk(const pst::Tree &st, const ptu::tree_meta &tm,
                     pt::idx_t m, pt::idx_t n, pt::idx_t ii_v_idx,
                     pt::idx_t ji_v_idx) {

    const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

    // std::cerr << pv_cmp::format("{}\n", fn_name);

    //auto [m, n] = mn;
    const std::vector<pt::idx_t> &height = tm.depth;
    // const std::vector<pt::idx_t> &lo = tm.lo;

    if (n == pc::INVALID_IDX || m == pc::INVALID_IDX || height[m] > height[n]) {
      return pc::INVALID_IDX; // invalid
    }

    pt::idx_t res = override_ji_trunk(st, tm, m, n, ii_v_idx, ji_v_idx);
    if (res != pc::INVALID_IDX) {
      return res;
    }

    // v_idx, height
    std::vector<std::pair<pt::idx_t, pt::idx_t>> x;

    for (pt::idx_t be_idx : st.get_obe_idxs(ji_v_idx)) {
      if (st.get_be(be_idx).type() != pst::be_type_e::back_edge) {
        continue;
      }

      pt::idx_t tgt_v_idx = st.get_be(be_idx).get_tgt();

      if (height[tgt_v_idx] < height[n]) {
        continue;
      }

      x.push_back({tgt_v_idx, height[tgt_v_idx]});
    }

    // std::cerr << fn_name << " x size: " << x.size() << "\n";

    // sort by height in ascending order
    std::sort(x.begin(), x.end(),
              [](const std::pair<pt::idx_t, pt::idx_t> &a,
                 const std::pair<pt::idx_t, pt::idx_t> &b) {
                return a.second < b.second;
              });

    if (x.empty()){
      return pc::INVALID_IDX;
    }
    else {
      return x.front().first;
    }


    for (auto [tgt_v_idx, h] : x) {
      bool valid{true};

      //std::cerr << fn_name << " tgt: " << tgt_v_idx << "\n";

      for (pt::idx_t be_idx : tm.get_brackets(tgt_v_idx)) {

        if (st.get_be(be_idx).type() != pst::be_type_e::back_edge) {
          continue;
        }

        pt::idx_t br_tgt_v_idx = st.get_be(be_idx).get_tgt();

        // if in trunk
        if (height[br_tgt_v_idx] > height[ii_v_idx]) {
          valid = false;
          break;
        }
      }

      if (!valid) {
        continue;
      }

      return tgt_v_idx;
    }

    return pc::INVALID_IDX;
}

void ji_branches(const pst::Tree &st, const ptu::tree_meta &tm,
                   pt::idx_t ii_v_idx, pt::idx_t ji_v_idx,
                   std::vector<pt::idx_t> &bb, pt::idx_t n) {

    // condition (i & ii)
    if (st.get_children(ji_v_idx).size() < 2 || ji_v_idx == n) {
      return;
    }

    pt::idx_t main_be_idx{pc::INVALID_IDX};
    for (pt::idx_t be_idx : st.get_ibe_idxs(ii_v_idx)) {
      const pst::BackEdge &be = st.get_be(be_idx);
      if (be.type() == pst::be_type_e::back_edge && be.get_src() == ji_v_idx) {
        main_be_idx = be_idx;
        break;
      }
    }

    const std::vector<pt::idx_t> &depth = tm.depth;
    // const std::vector<pt::idx_t> &height = tm.D;

    // children of j_i who meet condition (ii)
    std::vector<pt::idx_t> cond_iiia;
    std::vector<pt::idx_t> cond_iiib;
    for (pt::idx_t c_v_idx : st.get_children(ji_v_idx)) {
      pt::idx_t count{};
      bool has_br_into_ji{false};
      std::vector<pt::idx_t> br = tm.get_brackets(c_v_idx);
      for (pt::idx_t be_idx : br) {
        const pst::BackEdge &be = st.get_be(be_idx);
        pt::idx_t tgt_v_idx = be.get_tgt();

        if (tgt_v_idx != ji_v_idx) {
          count++;
        } else {
          has_br_into_ji = true;
        }
      }

      if (count == 1) {
        if (has_br_into_ji) {
          cond_iiia.push_back(c_v_idx);
        } else {
          cond_iiib.push_back(c_v_idx);
        }
      }
    }

    // condition (iii) (a)
    for (pt::idx_t c_v_idx : cond_iiia) {
      pt::idx_t lowest{c_v_idx};
      std::vector<pt::idx_t> br = tm.get_brackets(c_v_idx);
      for (pt::idx_t be_idx : br) {
        const pst::BackEdge &be = st.get_be(be_idx);
        pt::idx_t tgt_v_idx = be.get_tgt();
        pt::idx_t src_v_idx = be.get_src();

        if (tgt_v_idx != ji_v_idx) {
          continue;
        }

        if (depth[src_v_idx] > depth[lowest]) {
          lowest = src_v_idx;
        }
      }

#ifdef DEBUG
    assert(lowest != c_v_idx);
#endif

    bb.push_back(lowest);
  }

  auto is_branching = [&st](pt::idx_t v_idx) -> bool {
    return st.get_child_count(v_idx) > 1;
  };

  // condition (iii) (b)
  for (pt::idx_t c_v_idx : cond_iiib) {
    std::vector<pt::idx_t> brackets = tm.get_brackets(c_v_idx);

#ifdef DEBUG
    assert(brackets.size() == 1);
#endif

    pt::idx_t be_idx = brackets.front(); // size should be 1
    const pst::BackEdge &be = st.get_be(be_idx);
    pt::idx_t src_v_idx = be.get_src();
    pt::idx_t tgt_v_idx = be.get_tgt();

    // i. find d
    bool is_cond_i {false};
    pt::idx_t d {src_v_idx};
    while(d != ji_v_idx) {
      if (is_branching(d)) {
        is_cond_i = true;
        break;
      }
      d = st.get_parent_v_idx(d);
    }

    // ii. ensure for alpha
    //bool is_cond_ii {true};
    std::vector<pt::idx_t> alpha_br = tm.get_brackets(tgt_v_idx);

    if (!alpha_br.empty() && main_be_idx != pc::INVALID_IDX){
      continue;
    }

    if (alpha_br.size() != 1){
      continue;
    }

    if(is_cond_i){
      bb.push_back(d);
    }
  }
}

void with_ji(const pst::Tree &st, const ptu::tree_meta &tm,
          std::vector<pvst::Concealed> &res, pt::idx_t m, pt::idx_t n,
          pt::idx_t ii_v_idx, pt::idx_t ji_v_idx, pt::idx_t ft_v_idx) {
  const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  pt::idx_t tb = ji_trunk(st, tm, m, n, ii_v_idx, ji_v_idx);
  if (tb != pc::INVALID_IDX) {
    pvst::Concealed sl = gen_zi_slubble(st, ji_v_idx, tb, pvst::cl_e::zi_trunk, ft_v_idx);
    res.push_back(sl);
  }

  std::vector<pt::idx_t> bb; // branch boundaries
  ji_branches(st, tm, ii_v_idx, ji_v_idx, bb, n);
  for (auto b : bb) {
    pvst::Concealed sl = gen_zi_slubble(st, ji_v_idx, b, pvst::cl_e::zi_branch, ft_v_idx);
    res.push_back(sl);
  }
}
} // namespace zi


namespace update_pvst {

[[nodiscard]] inline bool is_nestable(pvst::Tree &pvst, pt::idx_t v_idx ) noexcept {
  const pvst::VertexBase &v = pvst.get_vertex(v_idx);
  return pvst::to_clan(v.get_fam()).value() == pvst::vc_e::fl_like;
  // return (v.get_type() == pvst::vt_e::flubble ||
  //                                    v.get_type() == pvst::vt_e::tiny ||
  //                                    v.get_type() == pvst::vt_e::parallel);
}

/*
  =================
  in relation to a
  =================
 */
void nest_trunk_ai(const pst::Tree &st, pvst::Tree &vst,
                   const ptu::tree_meta &tm, pt::idx_t sl_st_idx,
                   pt::idx_t fl_v_idx, pt::idx_t sl_v_idx,
                   std::vector<pt::idx_t> &ch) {
  const std::vector<pt::idx_t> &depth = tm.depth;

  for (pt::idx_t c_v_idx : ch) {

    if (!is_nestable(vst, c_v_idx)) {
      continue;
    }

    const pvst::Flubble &c_v = static_cast<const pvst::Flubble &>(vst.get_vertex(c_v_idx));

    pt::idx_t c_ai_idx = c_v.get_ai();
    pt::idx_t c_zi_idx = c_v.get_zi();

    bool above_g = depth[sl_st_idx] > depth[c_zi_idx];
    // bool from_g_branch =
    //     st.get_vertex(zi).post_order() < st.get_vertex(c_zi_idx).pre_order() &&
    //     st.get_vertex(sl_st_idx).pre_order() < st.get_vertex(c_ai_idx).pre_order() &&
    //     st.get_vertex(sl_st_idx).post_order() > st.get_vertex(c_zi_idx).post_order();
    bool from_g_branch = is_desc(st, sl_st_idx, c_ai_idx);

    if (above_g || from_g_branch) {
      vst.del_edge(fl_v_idx, c_v_idx);
      vst.add_edge(sl_v_idx, c_v_idx);
    }
  }
}

/**
 * @brief the fl lies between ℓ and a bracket of it goes into a_i
 */
void nest_branch_ai(const pst::Tree &st, pvst::Tree &vst,
                   const ptu::tree_meta &tm, pt::idx_t sl_st_idx,
                   pt::idx_t fl_v_idx, pt::idx_t sl_v_idx, pt::idx_t ai,
                   std::vector<pt::idx_t> &ch) {
  const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  // TODO: use pvst::bounds_t

  // is a descendant of...
  // a for ancestor, d for descendant
  // auto is_desc = [&](pt::idx_t a, pt::idx_t d) -> bool {
  //   return st.get_vertex(a).pre_order() < st.get_vertex(d).pre_order() &&
  //          st.get_vertex(a).post_order() > st.get_vertex(d).post_order();
  // };

  // could replace with v_idx.hi = a_i
  // does the vertex have a bracket into a_i
  auto has_ai_br = [&](pt::idx_t v_idx) -> bool {
    const std::vector<pt::idx_t> &br = tm.get_brackets(v_idx);
    for (pt::idx_t br_idx : br) {
      if (st.get_be(br_idx).get_src() == ai) {
        return true;
      }
    }
    return false;
  };

  for (pt::idx_t c_v_idx : ch) {

    if (!is_nestable(vst, c_v_idx)) {
      continue;
    }

    const pvst::Flubble &c_v = static_cast<const pvst::Flubble &>(vst.get_vertex(c_v_idx));

    pt::idx_t c_ai_idx = c_v.get_ai();
    pt::idx_t c_zi_idx = c_v.get_zi();

    bool desc = is_desc(st, sl_st_idx, c_ai_idx);
    bool has_br = has_ai_br(c_zi_idx);

    if (desc && has_br) {
      vst.del_edge(fl_v_idx, c_v_idx);
      vst.add_edge(sl_v_idx, c_v_idx);
    }
  }
}

void add_conc_ai(const pst::Tree &st, pvst::Tree &vst, const ptu::tree_meta &tm,
                 pt::idx_t fl_v_idx, const pvst::Flubble &v,
                 const std::vector<pvst::Concealed> &ai_adj, bool is_leaf) {

  const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  for (auto &sl : ai_adj) {
    pt::idx_t sl_v_idx = vst.add_vertex(sl);
    pt::idx_t sl_st_idx = sl.get_sl_st_idx();
    vst.add_edge(fl_v_idx, sl_v_idx);

    pt::idx_t ai = v.get_ai();

    if (!is_leaf) {
      // TODO: why not const vector ref?
      std::vector<pt::idx_t> ch = vst.get_children(fl_v_idx);
      if (sl.get_sl_type() == pvst::cl_e::ai_trunk) {
        nest_trunk_ai(st, vst, tm, sl_st_idx, fl_v_idx, sl_v_idx, ch);
      }
      else if (sl.get_sl_type() == pvst::cl_e::ai_branch) {
        nest_branch_ai(st, vst, tm, sl_st_idx, fl_v_idx, sl_v_idx, ai, ch);
      }
      else {
        std::cerr << "sl type: unknown\n";
      }
    }
  }
}

/*
  =================
  in relation to z
  =================
 */

/**
 * @brief if the flubble lies between n and z_i
 */
void nest_trunk_zi(const pst::Tree &st, pvst::Tree &vst,
                   const fl_sls &slubbles, pt::idx_t fl_v_idx,
                   pt::idx_t sl_v_idx, pt::idx_t zi,
                   std::vector<pt::idx_t> &ch) {
  const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  for (pt::idx_t c_v_idx : ch) {

    if (!is_nestable(vst, c_v_idx)) {
      continue;
    }

    const pvst::Flubble &c_v =
        static_cast<const pvst::Flubble &>(vst.get_vertex(c_v_idx));

    pt::idx_t c_ai_idx = c_v.get_ai();
    pt::idx_t c_zi_idx = c_v.get_zi();

    // a flubble is in the trunk if it a descendant of n but not a descendant of z_i
    bool desc_a = is_desc(st, slubbles.n, c_ai_idx);
    bool desc_b = !is_desc(st, zi, c_zi_idx);

    if (desc_a && desc_b) {
      vst.del_edge(fl_v_idx, c_v_idx);
      vst.add_edge(c_v_idx, sl_v_idx);
    }
  }
}

void nest_branch_zi(const pst::Tree &st, pvst::Tree &vst,
                   pt::idx_t sl_st_idx, pt::idx_t fl_v_idx, pt::idx_t sl_v_idx,
                   std::vector<pt::idx_t> &ch) {
  const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  for (pt::idx_t c_v_idx : ch) {
    if (!is_nestable(vst, c_v_idx)) {
      continue;
    }

    const pvst::Flubble &c_v =
        static_cast<const pvst::Flubble &>(vst.get_vertex(c_v_idx));

    // the fl is an ancestor of the slubble boundary
    bool ansc = is_desc(st, c_v.get_ai(), sl_st_idx);
    if (ansc) {
      vst.del_edge(fl_v_idx, c_v_idx);
      vst.add_edge(c_v_idx, sl_v_idx);
    }
  }
}

void add_conc_zi(const pst::Tree &st, const fl_sls &slubbles,
                 pvst::Tree &vst,
                 pt::idx_t fl_v_idx, const pvst::Flubble &v,
                 const std::vector<pvst::Concealed> &zi_adj, bool is_leaf) {
  const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  for (auto &sl : zi_adj) {
    pt::idx_t sl_v_idx = vst.add_vertex(sl);
    pt::idx_t sl_st_idx = sl.get_sl_st_idx();
    vst.add_edge(fl_v_idx, sl_v_idx);

    //pt::idx_t ai = v.get_ai();
    pt::idx_t zi = v.get_zi();

    if (!is_leaf) {
      // TODO: why not const vector ref?
      std::vector<pt::idx_t> ch = vst.get_children(fl_v_idx);

      if (sl.get_sl_type() == pvst::cl_e::zi_trunk) {
        nest_trunk_zi(st, vst, slubbles, fl_v_idx, sl_v_idx, zi, ch);
      }
      else if (sl.get_sl_type() == pvst::cl_e::ai_branch) {
        nest_branch_zi(st, vst, sl_st_idx, fl_v_idx, sl_v_idx, ch);
      }
      else {
        std::cerr << "sl type: unknown\n";
      }

    }
  }
}

/*
  =================

  =================
 */

/**
 * @brief update the PVST to add concealed bubbles
*/
void add_concealed(const pst::Tree &st, pvst::Tree &vst,
                   const ptu::tree_meta &tm, const fl_sls &slubbles) {
  const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  const auto &[fl_v_idx, ai_adj, zi_adj, _, __] = slubbles;

  const pvst::Flubble &v = static_cast<const pvst::Flubble &>(vst.get_vertex(fl_v_idx));

  bool is_leaf = st.get_children(fl_v_idx).empty();

  add_conc_ai(st, vst, tm, fl_v_idx, v, ai_adj, is_leaf);
  add_conc_zi(st, slubbles, vst, fl_v_idx, v, zi_adj, is_leaf);
}
} // namespace update_pvst

void find_concealed(const pst::Tree &st, pvst::Tree &ft, const ptu::tree_meta &tm) {
  const std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  std::vector<fl_sls> all_slubbles;
  for (pt::idx_t ft_v_idx{}; ft_v_idx < ft.vtx_count(); ft_v_idx++) {

    pvst::VertexBase &pvst_v = ft.get_vertex_mut(ft_v_idx);

    if (pvst_v.get_fam() != pvst::vt_e::flubble) {
      continue;
    }

    pvst::Flubble &ft_v = static_cast<pvst::Flubble &>(pvst_v);

    pt::idx_t ai = ft_v.get_ai();
    pt::idx_t zi = ft_v.get_zi();

    auto [m, n] = get_mn(st, tm, ai, zi);

    ft_v.set_m(m);
    ft_v.set_n(n);

    if (!can_contain(st, tm, ai, zi, m, n)) {
      continue;
    }

    // pass these as arguments to the slubble functions
    fl_sls slubbles;
    slubbles.fl_v_idx = ft_v_idx;
    slubbles.m = m;
    slubbles.n = n;
    ai::with_ai(st, tm, slubbles.ii_adj, m, n, ai, zi, ft_v_idx);
    zi::with_ji(st, tm, slubbles.ji_adj, m, n, ai, zi, ft_v_idx);


    if (slubbles.size() > 0) {

      all_slubbles.push_back(slubbles);
    }
  }

  for (auto &sl : all_slubbles) {
    update_pvst::add_concealed(st, ft, tm, sl);
  }
}
} // namespace povu::concealed
