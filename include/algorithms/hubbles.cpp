#include "./hubbles.hpp"
#include <cstddef>
#include <vector>

namespace povu::hubbles {

std::map<pt::idx_t, std::pair<pt::idx_t, pt::idx_t>> foo(const pst::Tree &st) {

  std::map<pt::idx_t, std::pair<pt::idx_t, pt::idx_t>> g_id_to_idx_map;

  /* spanning tree */
  for (pt::idx_t v_idx{}; v_idx < st.vtx_count(); v_idx++) {
    const pst::Vertex &v = st.get_vertex(v_idx);

    if (v.type() == pgt::v_type_e::dummy) {
      continue;
    }

    pt::idx_t g_v_id = v.g_v_id();

    if (v.type() == pgt::v_type_e::l) {

      if (!g_id_to_idx_map.contains(g_v_id)) {
        g_id_to_idx_map[g_v_id] = {v_idx, 0};
      } else {
        g_id_to_idx_map[g_v_id].first = v_idx;
      }
    }
    else if (v.type() == pgt::v_type_e::r) {
      if (!g_id_to_idx_map.contains(g_v_id)) {
        g_id_to_idx_map[g_v_id] = {0, v_idx};
      } else {
        g_id_to_idx_map[g_v_id].second = v_idx;
      }
    }
  }

  return g_id_to_idx_map;
}

struct fl_in_tr {
  pt::idx_t p_idx;
  pt::idx_t c_idx;
};

bool is_trunk_vtx(pst::Tree &st, pt::idx_t si_v_idx, pt::idx_t ei_v_idx,
                  pt::idx_t query_v_idx) {
  if (st.get_vertex(si_v_idx).post_order() > st.get_vertex(query_v_idx).post_order() &&
      st.get_vertex(query_v_idx).post_order() > st.get_vertex(ei_v_idx).post_order()) {
    return true;
  }

  return false;
}

std::vector<fl_in_tr> find_in_st(const pst::Tree &st,
                                     const pvtr::Tree<pgt::flubble> &ft) {

  std::map<pt::idx_t, std::pair<pt::idx_t, pt::idx_t>> g_id_to_idx_map = foo(st);

  std::vector<fl_in_tr> bar;

  /* flubble tree */
  for (pt::idx_t v_idx{}; v_idx < ft.vtx_count(); ++v_idx) {
    //const pgt::flubble &f = ft.get_vertex(v_idx).get_data().value();

    const pvtr::Vertex<pgt::flubble_t> &v = ft.get_vertex(v_idx);

    if (!v.get_data().has_value()) {
      continue;
    }

    auto [s, e] = v.get_data().value();
    auto [start_id, start_or] = s;
    auto [end_id, end_or] = e;

    pt::idx_t x = (start_or == pgt::or_e::forward)
        ? g_id_to_idx_map[start_id].second
        : g_id_to_idx_map[start_id].first;
    pt::idx_t y = (end_or == pgt::or_e::forward)
        ? g_id_to_idx_map[end_id].first
        : g_id_to_idx_map[end_id].second;

    if (st.get_vertex(x).dfs_num() < st.get_vertex(y).dfs_num()) {
      bar.push_back({x, y});
    }
    else {
      bar.push_back({y, x});
    }
  }

  return bar;
}

std::map<pt::idx_t, std::vector<std::pair<pt::idx_t, pt::idx_t>>>
annotate_branches(pst::Tree &st) {

  std::unordered_set<pt::idx_t> explored;

  std::stack<std::size_t> s;
  s.push( st.get_root_idx() );

  std::size_t prt_v_idx { st.get_root_idx() };
  std::size_t start { st.get_root_idx() };

  std::map<pt::idx_t, std::vector<std::pair<pt::idx_t, pt::idx_t>>> branch_map;
  std::map<pt::idx_t,  pt::idx_t> branch_map2;

  branch_map2[st.get_root_idx()] = st.vtx_count() - 1;

  while(!s.empty()) {
    pt::idx_t v_idx = s.top();

    if (explored.contains(v_idx)) {
      s.pop();
      continue;
    }

    const pst::Vertex &v = st.get_vertex(v_idx);

    if (v.is_leaf()) { // 0 children


      branch_map[prt_v_idx].push_back({start, v_idx});

      branch_map2[prt_v_idx] = v_idx;

      explored.insert(v_idx);
      s.pop();
      start = pc::INVALID_IDX;
    }
    else if (st.get_child_edges(v_idx).size() == 1) { // 1 child
      std::size_t c_v_idx = *st.get_children(v_idx).begin();
      s.push(c_v_idx);
      explored.insert(v_idx);
    }
    else if (st.get_children(v_idx).size() > 1) { // 2 or more, a branching path

      if (start !=pc::INVALID_IDX) {
        branch_map[prt_v_idx].push_back({start, v_idx});
      }

      prt_v_idx = v_idx;
      bool is_exp = true;
      for(auto c_v_idx : st.get_children(v_idx)) {
        if (!explored.contains(c_v_idx)) {
          start = c_v_idx;
          s.push(c_v_idx);
          is_exp = false;
          break;
        }
      }

      if (is_exp) {
        explored.insert(v_idx);
        start = pc::INVALID_IDX;
      }

    }
  }

  // for (auto [k, v] : branch_map) {
  //   std::cerr << std::format("branch: {} : ", k);
  //   for (auto [s, e] : v) {
  //     std::cerr << std::format("{} -> {}, ", s, e);
  //   }
  //   std::cerr << "\n";
  // }

  // for (auto [k, v] : branch_map2) {
  //   std::cerr << std::format("branch2: {} : ", k);
  //   std::cerr << std::format("{}\n", v);
  // }

  return branch_map;
}

struct tree_meta {
  std::vector<pt::idx_t> E;
  std::vector<pt::idx_t> D;
  //std::map<pt::idx_t, std::vector<std::pair<pt::idx_t, pt::idx_t>>> branch_map;
  std::vector<pt::idx_t> first; // idx is v_idx value is the first time it is seen in E
  std::vector<pt::idx_t> lo;

  void print() {

    // print lo
    std::cerr << "lo: \n";
    for (std::size_t v_idx = 0; v_idx < this->lo.size(); ++v_idx) {
      std::cerr << std::format("({}, {}), ", v_idx, this->lo[v_idx]);
    }
    std::cerr << "\n\n";

    // print E
    std::cerr << "E: \n";
    for (auto v_idx : this->E) {
      std::cerr << std::format("{} ", v_idx);
    }
    std::cerr << "\n\n";

    // print D
    std::cerr << "D: \n";
    for (std::size_t v_idx = 0; v_idx < this->D.size(); ++v_idx) {
      std::cerr << std::format("({}, {}), ", v_idx, this->D[v_idx]);
    }
    std::cerr << "\n\n";

    // print first
    std::cerr << "first (idx is v_idx value is the first time it is seen in E): \n";
    for (pt::idx_t v_idx = 0; v_idx < this->first.size(); ++v_idx) {
      std::cerr << std::format("({}, {}), ", v_idx, this->first[v_idx]);
    }
    std::cerr << "\n\n";
  }
};

/**
   *
   *
   *the max depth of a vertex reached by a backedge that starts below a
   * given backedge and ends at a vertex above this vertex
   */
void compute_lo(pst::Tree &st, tree_meta &tm) {

  std::vector<pt::idx_t> &lo = tm.lo;
  lo.reserve(st.vtx_count());
  for (pt::idx_t i = 0; i < st.vtx_count(); ++i) {
    lo.push_back(pc::INVALID_IDX);
  }

  for (pt::idx_t v_idx {st.vtx_count()} ; v_idx-- > 0 ; ) {

    pt::idx_t min_tgt{0};
    std::set<std::size_t> be_idxs = st.get_obe_idxs(v_idx);
    for (auto be_idx : be_idxs) {
      if (st.get_backedge(be_idx).type() != pst::be_type_e::back_edge) {
        // filter out special types of backedges
        continue;
      }

      pt::idx_t tgt_v_idx = st.get_backedge(be_idx).get_tgt();
      pt::idx_t tgt_depth = tm.D[tgt_v_idx];
      if (tgt_depth < tm.D[v_idx] && tgt_depth > min_tgt) {
        min_tgt = tgt_depth;
      }
    }

    pt::idx_t lo_child {0};
    for (pt::idx_t c_v_idx : st.get_children(v_idx)) {
      if (lo[c_v_idx] < tm.D[v_idx] && lo[c_v_idx] > lo_child) {
        lo_child = lo[c_v_idx];
      }
    }

    if (min_tgt > lo_child) {
      lo[v_idx] = min_tgt;
    }
    else {
      lo[v_idx] = lo_child;
    }
  }
}

void euler_tour(pst::Tree &st, tree_meta &tm) {

  std::unordered_set<pt::idx_t> explored;

  std::stack<std::size_t> s;
  s.push( st.get_root_idx() );

  std::vector<pt::idx_t> &E = tm.E;
  //std::vector<pt::idx_t> &D = tm.D;

  std::vector<pt::idx_t> D;
  for (pt::idx_t i = 0; i < st.vtx_count(); ++i) {
    D.push_back(0);
  }
  D[st.get_root_idx()] = 0;

  while(!s.empty()) {
    pt::idx_t v_idx = s.top();

    E.push_back(v_idx);

    if (explored.contains(v_idx)) {
      s.pop();
      continue;
    }

    const pst::Vertex &v = st.get_vertex(v_idx);

    if (!v.is_root()) {
      D[v_idx] = D[st.get_parent(v_idx)] + 1;
    }

    if (v.is_leaf()) { // 0 children
      explored.insert(v_idx);
      s.pop();
    }
    else if (st.get_child_edges(v_idx).size() == 1) { // 1 child
      std::size_t c_v_idx = *st.get_children(v_idx).begin();
      s.push(c_v_idx);
      explored.insert(v_idx);
    }
    else if (st.get_children(v_idx).size() > 1) { // 2 or more, a branching path
      bool is_exp = true;
      for(auto c_v_idx : st.get_children(v_idx)) {
        if (!explored.contains(c_v_idx)) {
          //start = c_v_idx;
          s.push(c_v_idx);
          is_exp = false;
          break;
        }
      }

      if (is_exp) {
        explored.insert(v_idx);
      }
    }
  }

  std::vector<pt::idx_t> &D_ = tm.D;
  D_.reserve(E.size());
  for (pt::idx_t i = 0; i < E.size(); ++i) {
    pt::idx_t v_idx = E[i];
    D_.push_back(D[v_idx]);
  }

  // first time we see a vertex in E
  std::vector<pt::idx_t> &first = tm.first;
  std::unordered_set<pt::idx_t> seen;

  for (pt::idx_t i = 0; i < E.size(); ++i) {
    pt::idx_t v_idx = E[i];
    if (!seen.contains(v_idx)) {
      seen.insert(v_idx);
      first.push_back(i);
    }
  }

  return ;
}



/**
 * put together backedges from the same branch of E_i
 */
std::map<pt::idx_t, std::vector<pt::idx_t>> cluster_be(pst::Tree &st,
                                                       pt::idx_t si_v_idx,
                                                       pt::idx_t ei_v_idx,
                                                       const std::vector<pt::idx_t> &branch_be_srcs) {

  std::map<pt::idx_t, std::vector<pt::idx_t>> branch_map;
  // TODO: [B] do in linear time
  for (auto ch_idx : st.get_children(ei_v_idx)) {
    const pst::Vertex& c_vtx = st.get_vertex(ch_idx);
    pt::idx_t pre_order = c_vtx.pre_order();
    pt::idx_t post_order = c_vtx.post_order();

    for (auto src_idx : branch_be_srcs) {
      const pst::Vertex& src_vtx = st.get_vertex(src_idx);
      pt::idx_t src_pre_order = src_vtx.pre_order();
      pt::idx_t src_post_order = src_vtx.post_order();

      if (src_pre_order >= pre_order && src_post_order <= post_order) {
        branch_map[ch_idx].push_back(src_idx);
      }
    }
  }

  return branch_map;
}

std::map<pt::idx_t, std::stack<pt::idx_t>>
compute_ancestors(pst::Tree &st, std::vector<pt::idx_t> &v, pt::idx_t limit) {

  /* find ancestors for the set of vertices and stop at limit */
  std::map<pt::idx_t, std::stack<pt::idx_t>> ancestor_map;

  //std::cerr << "limit " <<limit << "\n";

  //auto traverse_to_limit = [](){};

  for (auto src_v_idx : v) {
    pt::idx_t curr_v_idx = src_v_idx;

    /* traverse up the tree */
    while (true) {
      std::size_t p_v_idx = st.get_parent_v_idx(curr_v_idx);

      //std::cerr << p_v_idx <<", ";

      if (p_v_idx == limit) {
        //std::cerr << "\n";
        break;
      }

      //std::cerr << "pushing " << p_v_idx << ", ";

      ancestor_map[src_v_idx].push(p_v_idx);
      curr_v_idx = p_v_idx;
    }
  }

  return ancestor_map;
}



pt::idx_t find_lca(pst::Tree &st, const tree_meta &tm, std::vector<pt::idx_t> &vtxs, pt::idx_t ei_v_idx) {


  const std::vector<pt::idx_t> &first = tm.first;

  // TODO: [B] replace with O(1) op
  auto rmq = [&](pt::idx_t L, pt::idx_t R) {
    const auto &D = tm.D;
    pt::idx_t minDepth = std::numeric_limits<pt::idx_t>::max();
    pt::idx_t minPos = L;
    for (pt::idx_t i = L; i <= R; ++i) {
      if (D[i] < minDepth) {
        minDepth = D[i];
        minPos = i;
      }
    }
    return minPos;
  };

  // auto rmq = [&](pt::idx_t L, pt::idx_t R) {
  //   pt::idx_t min {pc::MAX_IDX};
  //   const std::vector<pt::idx_t> &D = tm.D;
  //   for (pt::idx_t i {L}; i <= R; i++){
  //     std::cerr << std::format("c {} {}, ", i, D[i]);
  //     if( D[i] < min) {
  //       min = D[i];
  //     }
  //   }

  //   return min;
  // };

  // 1) find the bounding interval [L…R]
  pt::idx_t L = {pc::MAX_IDX};
  pt::idx_t R = 0;
  for (pt::idx_t v : vtxs) {
    auto f = tm.first.at(v);
    L = std::min(L, f);
    R = std::max(R, f);
  }

  //pt::idx_t u = vtxs[0];
  //pt::idx_t v = vtxs[1];

  //std::cerr << "u: " << u << ", v: " << v << "\n";

  //pt::idx_t L = first.at(u);
  //pt::idx_t R = first.at(v);

  //std::cerr << "L: " << L << ", R: " << R << "\n";

  //if (L > R) {
  //   std::swap(L, R);
  //}

  pt::idx_t m = rmq(L, R);

  //std::cerr << "m: " << m << "\n";

  pt::idx_t lca = tm.E[m];

  return lca;
}



/**
 * from the branches
*/
std::vector<pt::idx_t> case_one(pst::Tree &st, const tree_meta &tm,
                                pt::idx_t si_v_idx, pt::idx_t ei_v_idx,
                                const std::vector<pt::idx_t> &branch_be_srcs) {


  std::vector<pt::idx_t> ends;

  /* has ancestor into trunk? */
  auto filter_be_to_trunk = [&](pt::idx_t lca_v_idx) -> bool {
    pt::idx_t curr_v_idx = lca_v_idx;

    //if (lca_v_idx == 7652) {}

    while (true) {

      // check if any of the obes of the curr vtx are going to the trunk
      std::set<std::size_t> obe_tgt_v_idxs = st.get_obe(curr_v_idx);
      for (auto tgt_v_idx : obe_tgt_v_idxs) {
        //pst::BackEdge be = st.get_backedge(obe_idx);

        //pt::idx_t tgt_v_idx = be.get_tgt();

        if (is_trunk_vtx(st, si_v_idx, ei_v_idx, tgt_v_idx)) {
          return true;
        }
      }

      // go up the tree
      if (st.get_parent_v_idx(curr_v_idx) == ei_v_idx) {
        break;
      } else {
        curr_v_idx = st.get_parent_v_idx(curr_v_idx);
      }
    }

    return false;
  };

  //std::cerr << "lca_v_idx: " << lca_v_idx << "\n";
  // cluster backedges from the same branch
  std::map<pt::idx_t, std::vector<pt::idx_t>> src_map =
    cluster_be(st, si_v_idx, ei_v_idx, branch_be_srcs);

  for (auto [k,v]: src_map) {
    if (v.size() > 1) {

      pt::idx_t lca_v_idx = find_lca(st, tm, v, ei_v_idx);

      if (!filter_be_to_trunk(lca_v_idx)) {
        ends.push_back(st.get_vertex(lca_v_idx).g_v_id());
      }

    }
  }

  return ends;
}

/**
 * from the trunk
*/
std::vector<pt::idx_t> case_two(pst::Tree &st,
                                pt::idx_t si_v_idx,
                                pt::idx_t ei_v_idx,
                                const std::vector<pt::idx_t>& trunk_be_srcs) {
  std::vector<pt::idx_t> ends;
  pt::idx_t boundary_vtx {};

  for (auto src_v_idx: trunk_be_srcs) {
    pt::idx_t curr_v_idx = src_v_idx;

    const pst::Vertex &v = st.get_vertex(curr_v_idx);

    if (v.is_leaf()) {
      // go up to branching path
      while (true) {
        std::size_t p_v_idx = st.get_parent_v_idx(curr_v_idx);

        const pst::Vertex &pv= st.get_vertex(p_v_idx);
        if (pv.get_obe().size() > 0) {
          boundary_vtx = pc::INVALID_IDX;
          break;
        }

        /* found a branching path */
        if (st.get_child_edges(p_v_idx).size() > 1) {
          boundary_vtx = p_v_idx;
          break;
        }
        curr_v_idx = p_v_idx;
      }
    }
    else {
      boundary_vtx = curr_v_idx;
    }

    if (boundary_vtx != pc::INVALID_IDX) {
      ends.push_back(st.get_vertex(boundary_vtx).g_v_id());
    }

  }

  return ends;
}

void filter_branches(pst::Tree &st, pt::idx_t si_v_idx, pt::idx_t ei_v_idx) {

  std::map<pt::idx_t, pt::idx_t> branches;

  for (pt::idx_t c_v_idx : st.get_children(ei_v_idx)) {
    branches[c_v_idx] = pc::INVALID_IDX;
  }

}

std::pair<std::vector<pt::idx_t>, std::vector<pt::idx_t>>
split_back_edges(pst::Tree &st, const tree_meta &tm, pt::idx_t si_v_idx,
                 pt::idx_t ei_v_idx) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  std::vector<pt::idx_t> trunk_be_srcs;
  std::vector<pt::idx_t> branch_be_srcs;

  std::set<std::size_t> ibe_idxs = st.get_ibe_idxs(si_v_idx);
  for (auto ibe_idx : ibe_idxs) {
    pst::BackEdge be = st.get_backedge(ibe_idx);

    if (be.type() != pst::be_type_e::back_edge) {
      continue;
    }

    pt::idx_t src_v_idx = be.get_src();

    if (st.get_vertex(src_v_idx).post_order() > st.get_vertex(ei_v_idx).post_order()) {
      // starts from a trunk to si
      trunk_be_srcs.push_back(src_v_idx);
    }
    else {
      branch_be_srcs.push_back(src_v_idx);
    }
  }

  std::vector<pt::idx_t> x;

  std::set<std::size_t> child_v_idxs = st.get_children(ei_v_idx);
  for (std::size_t c_v_idx: child_v_idxs){

    if (st.get_vertex(c_v_idx).hi() == si_v_idx && tm.lo[c_v_idx] == si_v_idx) {
      for(pt::idx_t src_v_idx : branch_be_srcs) {
        if (st.get_vertex(c_v_idx).pre_order() < st.get_vertex(src_v_idx).pre_order() &&
            st.get_vertex(c_v_idx).post_order() > st.get_vertex(src_v_idx).post_order()) {
          x.push_back(src_v_idx);
        }
      }
    }
  }


  std::pair<std::vector<pt::idx_t>, std::vector<pt::idx_t>> be_srcs{ trunk_be_srcs, x };
  return be_srcs;
}

std::vector<pt::idx_t> with_si(pst::Tree &st, const tree_meta &tm, pt::idx_t si_v_idx, pt::idx_t ei_v_idx) {

  pt::idx_t si_v_id = st.get_vertex(si_v_idx).g_v_id();

  /* split the set of back edges into trunk and branch backedges */
  auto [trunk_be_srcs, brach_be_srcs] = split_back_edges(st, tm, si_v_idx, ei_v_idx);

  std::vector<pt::idx_t> ends_one = case_one(st, tm, si_v_idx, ei_v_idx, brach_be_srcs);

  for (auto e : ends_one) {

    std::cerr << std::format("branch boundary: {} {}\n", si_v_id, e);
  }


  std::vector<pt::idx_t> ends = case_two(st, si_v_idx, ei_v_idx, trunk_be_srcs);
  for (auto e: ends) {
    std::cerr << std::format("trunk boundary: {} {}\n", si_v_id, e);
  }

  ends.insert(ends.end(), ends_one.begin(), ends_one.end());

  return ends;
}

// debug fn
void foo(pst::Tree &st, const pvtr::Tree<pgt::flubble> &ft) {
  std::vector<fl_in_tr> fl = find_in_st(st, ft);

  for (auto [si, ei] : fl) {
    if (st.get_vertex(si).g_v_id() >= 1867 && st.get_vertex(ei).g_v_id() <= 1873) {

      // print src of ibe into si
      //const pst::Vertex &v =st.get_ibe_idxs(si);

      std::cerr <<std::format("si: {} ei: {}\n", si, ei);

      std::set<std::size_t> ibe_idxs = st.get_ibe_idxs(si);

      for (auto ibe_idx : ibe_idxs) {
        pst::BackEdge be = st.get_backedge(ibe_idx);

        if (be.type() != pst::be_type_e::back_edge) {
          continue;
        }

        pt::idx_t src_v_idx = be.get_src();

        std::cerr <<std::format("src: {} {} \n", src_v_idx, st.get_vertex(src_v_idx).g_v_id());

        // for each src print parent until ei
        pt::idx_t curr_v_idx = src_v_idx;
        while (true) {
          std::size_t p_v_idx = st.get_parent_v_idx(curr_v_idx);


          std::cerr <<std::format("parent: {} {} \n", p_v_idx, st.get_vertex(p_v_idx).g_v_id());

          if (p_v_idx == ei) {
            break;
          }

          curr_v_idx = p_v_idx;
        }
      }

      std::cerr << "ei prts \n";


      // print parent from ei to si
      pt::idx_t curr_v_idx = ei;
      while (true) {
        std::size_t p_v_idx = st.get_parent_v_idx(curr_v_idx);

        std::cerr <<std::format("parent: {} {} \n", p_v_idx, st.get_vertex(p_v_idx).g_v_id());

        if (p_v_idx == si) {
          break;
        }

        curr_v_idx = p_v_idx;
      }
    }
  }

}

void find_hubbles(pst::Tree &st, const pvtr::Tree<pgt::flubble> &ft) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  // debug fn
  foo(st, ft);

  std::vector<fl_in_tr> bar = find_in_st(st, ft);
  tree_meta tm;
  euler_tour(st, tm);
  compute_lo(st, tm);
  //tm.print();
  /* a for ancestor, d for descendant */
  for (auto [si, ei] : bar) {
    //std::cerr << std::format("x: {} y: {}\n", si, ei);
    //pt::idx_t si_v_id = st.get_vertex(si).g_v_id();
    //pt::idx_t boundary = with_top(st, si, ei);
    for (auto _ : with_si(st, tm, si, ei)) {
      //std::cerr << std::format("boundary: {} {}\n", si_v_id, b);
    }
    //std::cerr << std::format("boundary: {} {}\n", si_v_id, boundary);
  }
}

} // namespace povu::hubbles
