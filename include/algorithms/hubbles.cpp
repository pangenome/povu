#include "./hubbles.hpp"
#include <cstddef>
#include <iostream>
#include <queue>
#include <unordered_map>
#include <utility>
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
  pt::idx_t ii_idx; // parent idx
  pt::idx_t ji_idx; // child idx
};

bool is_trunk_vtx(pst::Tree &st, pt::idx_t si_v_idx, pt::idx_t ei_v_idx,
                  pt::idx_t query_v_idx) {
  if (st.get_vertex(si_v_idx).post_order() > st.get_vertex(query_v_idx).post_order() &&
      st.get_vertex(query_v_idx).post_order() > st.get_vertex(ei_v_idx).post_order()) {
    return true;
  }

  return false;
}

std::vector<fl_in_tr> find_in_st(const pst::Tree &st, const pvtr::Tree<pgt::flubble> &ft) {

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

/*
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
*/

struct tree_meta {
  std::vector<pt::idx_t> E;
  std::vector<pt::idx_t> D;
  //std::map<pt::idx_t, std::vector<std::pair<pt::idx_t, pt::idx_t>>> branch_map;
  std::vector<pt::idx_t> first; // idx is v_idx value is the first time it is seen in E
  std::vector<pt::idx_t> lo;

  std::map<pt::idx_t, pt::idx_t> pre; // idx is the pre-order the value is the v_idx
  std::map<pt::idx_t, pt::idx_t> post; // idx is the post-order the value is the v_idx

  // Gather all backedges
  std::vector<pt::idx_t> B;

  // prefix sum of the number of backedges
  std::vector<pt::idx_t> off;

  // a flat list of backedges
  std::vector<pt::idx_t> BE;

  std::vector<pt::idx_t> get_brackets(pt::idx_t v_idx) const {
    std::vector<pt::idx_t> brackets;
    pt::idx_t start = off[v_idx];
    pt::idx_t end = off[v_idx + 1];

    for (pt::idx_t i{start}; i < end; i++) {
      pt::idx_t be_idx = BE[i];
      brackets.push_back(be_idx);
    }

    return brackets;
  }

  std::vector<pt::idx_t> depth;

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

    // print depth
    std::cerr << "depth: \n";
    for (pt::idx_t v_idx = 0; v_idx < this->depth.size(); ++v_idx) {
      std::cerr << std::format("({}, {}), ", v_idx, this->depth[v_idx]);
    }
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

void compute_bracket_props(pst::Tree &st) {
  for (pt::idx_t v_idx{st.vtx_count()}; v_idx-- > 0;) {
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



pt::idx_t find_lca(const tree_meta &tm, std::vector<pt::idx_t> &vtxs) {


  //const std::vector<pt::idx_t> &first = tm.first;

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

/*
  -----------------
  Handle S_i
  -----------------
*/




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

struct mn_t {
  pt::idx_t m;
  pt::idx_t n;
};

pt::idx_t compute_m(pst::Tree &st, const tree_meta &tm, pt::idx_t ii_v_idx, pt::idx_t ji_v_idx) {
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
    if (st.get_backedge(be_idx).type() != pst::be_type_e::back_edge) {
      continue;
    }

    pt::idx_t src_v_idx = st.get_backedge(be_idx).get_src();

    if (height[src_v_idx] >= height[ji_v_idx]) {
      continue;
    }

    if (height[src_v_idx] > max_height.height) {
      // TODO: [c] remove type cast
      max_height = {(pt::idx_t)src_v_idx, height[src_v_idx]};
    }
  }

  return  max_height.v_idx == pc::INVALID_IDX ? ii_v_idx : max_height.v_idx;
}

pt::idx_t compute_n(pst::Tree &st, const tree_meta &tm, pt::idx_t ii_v_idx, pt::idx_t ji_v_idx) {

  pt::idx_t count{0};
  // count normal OBE
  for (pt::idx_t be_idx : st.get_obe_idxs(ji_v_idx)) {
    if (st.get_backedge(be_idx).type() == pst::be_type_e::back_edge) {
      count++;
    }
  }

  if (count == 0) {
    return ji_v_idx;
  }

  //std::cerr << "OBE cnt " << st.get_obe_tgt_v_idxs(ji_v_idx).size() << "\n";

  const std::vector<pt::idx_t> &height = tm.depth;

  struct v_idx_height_t {
    pt::idx_t v_idx;
    pt::idx_t height;
  };

  v_idx_height_t min_height{pc::INVALID_IDX, pc::MAX_IDX};
  for (auto tgt_v_idx : st.get_obe_tgt_v_idxs(ji_v_idx)) {
    if (ii_v_idx == tgt_v_idx) {
      continue;
    }

    if (height[tgt_v_idx] < min_height.height) {
      // TODO: [c] remove type cast
      min_height = {(pt::idx_t)tgt_v_idx, height[tgt_v_idx]};
    }
  }

  return  min_height.v_idx;
}

/**
 */
mn_t get_mn(pst::Tree &st, const tree_meta &tm, pt::idx_t ii_v_idx,
            pt::idx_t ji_v_idx) {
  // when m is invalid (i) of ii trunk be is false

  pt::idx_t m = compute_m(st, tm, ii_v_idx, ji_v_idx);
  pt::idx_t n = compute_n(st, tm, ii_v_idx, ji_v_idx);

  return {m, n};
}

/**
get Ii trunk back edges
 */
pt::idx_t ii_trunk(pst::Tree &st, const tree_meta &tm, const mn_t &mn,
                   pt::idx_t ii_v_idx, pt::idx_t ji_v_idx) {

  auto [m, n] = mn;
  std::vector<pt::idx_t> height = tm.depth;

  if (n == pc::INVALID_IDX || m == pc::INVALID_IDX || height[m] > height[n]) {
    return pc::INVALID_IDX; // invalid
  }

  // v_idx, height
  std::vector<std::pair<pt::idx_t, pt::idx_t>> x;

  for (auto be_src_v_idx : st.get_ibe_src_v_idxs(ii_v_idx)) {
    if (height[be_src_v_idx] > height[m]) {
      continue;
    }

    x.push_back({be_src_v_idx, height[be_src_v_idx]});
  }

  // sort by height in descending order
  std::sort(x.begin(), x.end(),
            [](const std::pair<pt::idx_t, pt::idx_t> &a,
               const std::pair<pt::idx_t, pt::idx_t> &b) {
              return a.second > b.second;
            });

  // loop from start to end of x and get the first one whose bracket source is
  // not from below m and is not ji already handled in x
  for (auto [be_src_v_idx, h] : x) {
    bool invalid {false};
    for (auto be_idx : tm.get_brackets(be_src_v_idx)) {
      const pst::BackEdge &be = st.get_backedge(be_idx);
      pt::idx_t src_v_idx = be.get_src();

      if (height[src_v_idx] > height[m]) {
        invalid = true;
        break;
      }
    }

    if (invalid) {
      continue;
    }

    std::vector<pt::idx_t> s = { be_src_v_idx, ji_v_idx };
    return find_lca(tm, s);
    //return be_src_v_idx;
    //return st.get_vertex(be_src_v_idx).g_v_id();
  }

  return pc::INVALID_IDX;
}

/** get Ii brach backedges */
void ii_branches(pst::Tree &st, const tree_meta &tm,
                 pt::idx_t ii_v_idx, pt::idx_t ji_v_idx,
                 std::vector<pt::idx_t> &bb) {

  if (st.get_children(ji_v_idx).size() < 2) {
    return;
  }

  // std::map<pt::idx_t, pt::idx_t> branches;

  const std::vector<pt::idx_t> &height = tm.D;
  const std::vector<pt::idx_t> &lo = tm.lo;

  // children who meet condition (i) and (ii)
  std::vector<pt::idx_t> x;
  for (pt::idx_t c_v_idx : st.get_children(ji_v_idx)) {
    if (st.get_vertex(c_v_idx).hi() == lo[c_v_idx] && st.get_vertex(c_v_idx).hi() == ii_v_idx) {
      x.push_back(c_v_idx);
    }
  }

  for (auto c_v_idx : x) {

    // std::cerr << "child" << c_v_idx << "\n";

    std::vector<pt::idx_t> brackets = tm.get_brackets(c_v_idx);

    // case (ii) (a)
    if (brackets.size() < 2) {
      continue;
    }

    // std::cerr << "bracket count: "<< brackets.size() << "\n";

    std::vector<pt::idx_t> br_srcs;
    for (auto be_idx : brackets) {
      const pst::BackEdge &be = st.get_backedge(be_idx);
      pt::idx_t src_v_idx = be.get_src();
      br_srcs.push_back(src_v_idx);
    }

    pt::idx_t d = find_lca(tm, br_srcs);
    if (tm.get_brackets(d).size() > 0) {
      bb.push_back(st.get_vertex(d).g_v_id());
    }
    else {
      // get the br srcs with min depth
      // v_idx, height
      // pt::idx_t min;
      std::pair<pt::idx_t, pt::idx_t> min_depth {pc::INVALID_IDX, pc::MAX_IDX};
      for (pt::idx_t be_idx : tm.get_brackets(d)) {
        const pst::BackEdge &be = st.get_backedge(be_idx);
        pt::idx_t br_src = be.get_src();
        if (height[br_src] < min_depth.second) {
          min_depth = {br_src, height[br_src]};
        }
      }

      if (min_depth.first == pc::INVALID_IDX) {
        std::cerr << "invalid min for lca: " << d << "of: ";
        for (auto x : br_srcs){
          std::cerr << x << ", ";
        }
        std::cerr << "\n";
      }
      else {
        //std::cerr << "min " << min_depth.first << "\n";
        bb.push_back(st.get_vertex(min_depth.first).g_v_id());
      }

    }
  }

}

std::vector<pt::idx_t> with_ii(pst::Tree &st, const tree_meta &tm,
                               const mn_t &mn, pt::idx_t ii_v_idx,
                               pt::idx_t ji_v_idx) {


  pt::idx_t tb = ii_trunk(st, tm, mn, ii_v_idx, ji_v_idx);
  if (tb != pc::INVALID_IDX) {
    std::cerr << std::format("ii trunk boundary: {} {}\n",
                             st.get_vertex(ii_v_idx).g_v_id(),
                             st.get_vertex(tb).g_v_id());
  }
  //std::cerr << std::format("ii trunk boundary: {} {}\n",
  //                         st.get_vertex(ii_v_idx).g_v_id(),
  //                         st.get_vertex(tb).g_v_id());

  std::vector<pt::idx_t> bb; // branch boundaries
  ii_branches(st, tm, ii_v_idx, ji_v_idx, bb);
  for (auto b : bb) {
    std::cerr << std::format("ii branch boundary: {} {}\n",
                             st.get_vertex(ii_v_idx).g_v_id(),
                             b);
  }

  return bb;
}

/*
  -----------------
  Handle E_i
  -----------------
*/

pt::idx_t ji_trunk(pst::Tree &st, const tree_meta &tm, const mn_t &mn,
                   pt::idx_t ii_v_idx, pt::idx_t ji_v_idx) {

  auto [m, n] = mn;
  const std::vector<pt::idx_t> &height = tm.D;
  const std::vector<pt::idx_t> &lo = tm.lo;

  if (n == pc::INVALID_IDX || m == pc::INVALID_IDX || height[m] > height[n]) {
    return pc::INVALID_IDX; // invalid
  }

  // v_idx, height
  std::vector<std::pair<pt::idx_t, pt::idx_t>> x;

  for (pt::idx_t tgt_v_idx : st.get_obe_tgt_v_idxs(ji_v_idx)) {
    if (height[tgt_v_idx] < height[n]) {
      continue;
    }
    x.push_back({tgt_v_idx, height[tgt_v_idx]});
  }

  // sort by height in ascending order
  std::sort(x.begin(), x.end(),
            [](const std::pair<pt::idx_t, pt::idx_t> &a,
               const std::pair<pt::idx_t, pt::idx_t> &b) {
              return a.second < b.second;
            });

  for (auto [tgt_v_idx, h] : x) {
    if (lo[tgt_v_idx] >= height[tgt_v_idx]) {
      return tgt_v_idx;
    }
  }

  return pc::INVALID_IDX;
}

void ji_branches(pst::Tree &st, const tree_meta &tm, pt::idx_t ii_v_idx,
                 pt::idx_t ji_v_idx, std::vector<pt::idx_t> &bb) {

  // condition (i)
  if (st.get_children(ji_v_idx).size() < 2) {
    return;
  }

  const std::vector<pt::idx_t> &height = tm.D;

  // children of j_i who meet condition (ii)
  std::vector<pt::idx_t> x;
  // get children with only one child into the trunk
  for (pt::idx_t c_v_idx : st.get_children(ji_v_idx)) {
    bool valid {true};
    pt::idx_t count {0}; // number of trunk vertices
    for (pt::idx_t be_idx: tm.get_brackets(c_v_idx)) {
      const pst::BackEdge &be = st.get_backedge(be_idx);
      pt::idx_t tgt_v_idx = be.get_tgt();

      if (tgt_v_idx == ii_v_idx) {
        valid= false;
        break;
      }

      if (height[tgt_v_idx] < height[ji_v_idx]) {
        count++;
      }

      if (count > 1) {
        valid = false;
        break;
      }
    }

    if (valid) {
      x.push_back(c_v_idx);
    }
  }

  std::vector<pt::idx_t> y; // children who meet condition (iv) (b)

  // condition (iv) (a)
  for (pt::idx_t c_v_idx : x) {
    pt::idx_t max {};
    for (pt::idx_t be_idx : tm.get_brackets(c_v_idx)) {
      const pst::BackEdge &be = st.get_backedge(be_idx);
      pt::idx_t tgt_v_idx = be.get_tgt();
      pt::idx_t src_v_idx = be.get_src();
      if (tgt_v_idx == ji_v_idx) {
        if (height[src_v_idx] > height[max]) {
          max = src_v_idx;
        }
      }
    }

    if (max == 0) {
      y.push_back(c_v_idx);
    }
    else {
      // meet condition (iv) (a)
      std::cerr << "pushing back (iv) f(c) " << max << "\n";
      bb.push_back(max);
    }
  }

  // condition (iv) g(c)
  for (pt::idx_t c_v_idx : y) {
    std::vector<pt::idx_t> brackets = tm.get_brackets(c_v_idx);

    if (brackets.size() != 1) {
      //std::cerr << std::format("handle case brackets size: {} \n", brackets.size());
      continue;
    }

    //if (brackets.size() != 1) {

    //}


    std::cerr << "child" << c_v_idx << "\n";


    std::cerr << "a\n";
    pt::idx_t be_idx = brackets.front(); // size should be 1
    std::cerr << "b\n";
    const pst::BackEdge &be = st.get_backedge(be_idx);
    pt::idx_t src_v_idx = be.get_src();
    pt::idx_t tgt_v_idx = be.get_tgt();

    std::cerr << src_v_idx << " - - -> " << tgt_v_idx << "\n";

    pt::idx_t curr_v_idx = src_v_idx;
    while (st.get_parent_v_idx(curr_v_idx) != ji_v_idx) {
      std::size_t p_v_idx = st.get_parent_v_idx(curr_v_idx);

      if (st.get_child_edges(p_v_idx).size() > 1) {
        std::cerr << "pushing back (iv) g(c) " << p_v_idx << "\n";
        bb.push_back(p_v_idx);
        break;
      }

      curr_v_idx = p_v_idx;
    }

    // find lowest branching ancestor that is not a trunk vertex (including ji)

    // while (curr_v_idx != ji_v_idx) {
    //   std::size_t p_v_idx = st.get_parent_v_idx(curr_v_idx);

    //   if (st.get_child_edges(p_v_idx).size() > 1) {
    //     std::cerr << "pushing back (iv) g(c) " << p_v_idx << "\n";
    //     bb.push_back(p_v_idx);
    //     break;
    //   }

    //   curr_v_idx = p_v_idx;
    // }
  }
}

// with trunk
pt::idx_t ei_case_one(pst::Tree &st, const tree_meta &tm,
                                    pt::idx_t si_v_idx, pt::idx_t ei_v_idx,
                                   pt::idx_t min_depth_v_idx) {

  std::vector<pt::idx_t> ends;

  const std::vector<pt::idx_t> &depth = tm.D;

  // max depth
  struct max_d {
    pt::idx_t v_idx;
    pt::idx_t depth;
  };
  max_d min_depth{pc::INVALID_IDX, pc::INVALID_IDX};

  for (pt::idx_t be_idx : st.get_obe_idxs(ei_v_idx)) {
    const pst::BackEdge& be = st.get_backedge(be_idx);

    if (be.type() != pst::be_type_e::back_edge) {
      continue;
    }

    pt::idx_t tgt_v_idx = be.get_tgt();
    //std::cerr << std::format("tgt_v_idx: {} \n", tgt_v_idx);
    if (depth[tgt_v_idx] > depth[min_depth_v_idx] && depth[tgt_v_idx] < min_depth.depth) {
      min_depth = {tgt_v_idx, tm.D[tgt_v_idx]};
    }
  }

  // print ends
  if (min_depth.v_idx != pc::INVALID_IDX) {
    ends.push_back(st.get_vertex(min_depth.v_idx).g_v_id());

     std::cerr << "with E_i\n";
    // std::cerr << "si " << si_v_idx << " ei " << ei_v_idx << "\n";
    //std::cerr << min_depth_v_idx << " " << ei_v_idx << "\n";
    std::cerr << "v " << st.get_vertex(min_depth.v_idx).g_v_id() << " E_i " << st.get_vertex(ei_v_idx).g_v_id() << "\n";
    std::cerr << "------------\n";
  }


  return min_depth_v_idx;
}

std::vector<pt::idx_t> ei_case_two(pst::Tree &st, const tree_meta &tm,
                                   pt::idx_t si_v_idx, pt::idx_t ei_v_idx,
                                   pt::idx_t min_depth_v_idx, std::vector<pt::idx_t> &ends) {

  const std::vector<pt::idx_t> &depth = tm.D;
  const std::vector<pt::idx_t> &lo = tm.lo;

  std::vector<pt::idx_t> ch;

  // filter out children whose hi is not less than min depth
  for (pt::idx_t c_v_idx : st.get_children(ei_v_idx)) {
    pt::idx_t hi_v_idx = st.get_vertex(c_v_idx).hi();
    bool hi_is_below_min = depth[hi_v_idx] > depth[min_depth_v_idx];
    bool hi_is_above_ei = depth[hi_v_idx] < depth[ei_v_idx];

    if (hi_is_above_ei && hi_is_below_min ) {
      ch.push_back(c_v_idx);
    }
  }

  // get only children descendants have 1 obe and the tgt of that OBE lies between min depth and
  // E_i

  // only one bracket must go above E_i
  std::set<pt::idx_t> invalid_c;
  for (pt::idx_t c_v_idx : ch) {
    std::vector<pt::idx_t> brackets = tm.get_brackets(c_v_idx);
    //bool is_valid = true;
    pt::idx_t cnt {}; // no of brackets that go above E_i
    for (auto be_idx : brackets) {
      if(cnt > 1){
        invalid_c.insert(c_v_idx);
        break;
      }

      const pst::BackEdge& be = st.get_backedge(be_idx);


      pt::idx_t y = be.get_tgt();

      if (depth[y] < depth[ei_v_idx]) {
        cnt++;
      }
     }
  }

  // ch without be into E_i
  for  (auto c_v_idx : ch) {

    if(invalid_c.contains(c_v_idx)) {
      continue;
    }

    // TODO:[B] query the count
    std::vector<pt::idx_t> brackets = tm.get_brackets(c_v_idx);
    if (brackets.size() > 1) {
      continue;
    }

    pt::idx_t be_idx = brackets[0];
    const pst::BackEdge& be = st.get_backedge(be_idx);
    pt::idx_t x = be.get_src();
    pt::idx_t y = be.get_tgt();

    //std::cerr << std::format("src: {} tgt: {} \n", st.get_vertex(x).g_v_id(), st.get_vertex(y).g_v_id());

    if (!st.is_leaf(x)){
      continue;
    }

    // the only be we allow to bracket y between E_i and S_i is the E_i, S_i back edge
    bool y_is_invalid = false;
    for (auto be_idx : tm.get_brackets(y)){
      const pst::BackEdge& be = st.get_backedge(be_idx);
      pt::idx_t src_v_idx = be.get_src();
      pt::idx_t tgt_v_idx = be.get_tgt();

      bool a = src_v_idx == ei_v_idx && tgt_v_idx == si_v_idx; // allowed

      bool b = depth[src_v_idx] < depth[ei_v_idx] || depth[tgt_v_idx] > depth[si_v_idx]; //
      // not allowed

      if (!a || b) {
        y_is_invalid = true;
        break;
      }
    }

    if (y_is_invalid) {
      continue;
    }

    std::cerr << std::format("with E_i branch  ({}: {}, {})\n",
                             st.get_vertex(c_v_idx).g_v_id(),
                             st.get_vertex(y).g_v_id(),
                             st.get_vertex(ei_v_idx).g_v_id());
  }


  // ch with be into E_i
  // get the lowest src of a be into E_i
  for  (auto c_v_idx : ch) {

    if(invalid_c.contains(c_v_idx)) {
      continue;
    }

    bool lo_is_ei = lo[c_v_idx] == ei_v_idx; // we have a be into E_i

    if (!lo_is_ei) {
      continue;
    }

    pt::idx_t max_src {c_v_idx};

    // TODO:[B] query the count
    std::vector<pt::idx_t> brackets = tm.get_brackets(c_v_idx);
    for (auto be_idx : brackets) {
      const pst::BackEdge &be = st.get_backedge(be_idx);

      pt::idx_t x = be.get_src();
      pt::idx_t y = be.get_tgt();

      //std::cerr << std::format("br src: {} tgt: {} \n", st.get_vertex(x).g_v_id(), st.get_vertex(y).g_v_id());

      if (y != ei_v_idx) {
        continue;
      }



      if (x > max_src) {
        max_src = x;
      }

    }

    //pt::idx_t be_idx = brackets[0];
    //const pst::BackEdge& be = st.get_backedge(be_idx);

    //pt::idx_t x = be.get_src();

    // x has a descendant that whose backedge goes to E_i

    //std::cerr << "is leaf " << st.is_leaf(max_src) << "\n";

    pt::idx_t a;
    if (st.is_leaf(max_src)) {
      a = st.get_vertex(c_v_idx).hi();
    }
    else {
      a = max_src;
    }


    // note that the closest thing to a flubble in this case is the lowest src of a be into E_i

    ends.push_back(st.get_vertex(max_src).g_v_id());


    std::cerr << std::format("with E_i branch be into E_i ({} {}, {})\n",
                              st.get_vertex(c_v_idx).g_v_id(),
                             st.get_vertex(a).g_v_id(),
                             st.get_vertex(ei_v_idx).g_v_id());
    //std::cerr << "----\n";
  }

  return ends;
}

// get the max depth of a backedge from the trunk into S_i
pt::idx_t find_min_depth_vtx(pst::Tree &st, const tree_meta &tm, pt::idx_t si_v_idx, pt::idx_t ei_v_idx) {

  const std::vector<pt::idx_t> &depth = tm.D;

  // max depth
  struct max_d {
    pt::idx_t v_idx;
    pt::idx_t depth;
  };
  max_d max_depth{pc::INVALID_IDX, 0};

  for (pt::idx_t be_idx : st.get_ibe_idxs(si_v_idx)) {
    const pst::BackEdge& be = st.get_backedge(be_idx);

    if (be.type() != pst::be_type_e::back_edge) {
      continue;
    }

    pt::idx_t src_v_idx = be.get_src();
    // excludes backedge E_i -> S_i
    if (depth[src_v_idx] > max_depth.depth && depth[src_v_idx] < depth[ei_v_idx]) {
      max_depth = {src_v_idx, depth[src_v_idx]};
    }
  }

  // happens if there are no backedges or if the only backedge is E_i -> S_i
  if (max_depth.v_idx == pc::INVALID_IDX){
    max_depth.v_idx = si_v_idx;
  }

  return max_depth.v_idx;
}

std::vector<pt::idx_t> with_ji(pst::Tree &st, const tree_meta &tm,
                               const mn_t &mn, pt::idx_t ii_v_idx,
                               pt::idx_t ji_v_idx) {
  std::vector<pt::idx_t> ends;

  pt::idx_t tb = ji_trunk(st, tm, mn, ii_v_idx, ji_v_idx);
  if (tb != pc::INVALID_IDX) {
    ends.push_back(st.get_vertex(tb).g_v_id());
    std::cerr << std::format("ji trunk boundary: {} {}\n",
                             st.get_vertex(ji_v_idx).g_v_id(),
                             st.get_vertex(tb).g_v_id());
  }
  //std::cerr << std::format("ji trunk boundary: {} {}\n",
  //                         st.get_vertex(tb).g_v_id(),
  //                         st.get_vertex(ji_v_idx).g_v_id());

  std::vector<pt::idx_t> bb; // branch boundaries
  ji_branches(st, tm, ii_v_idx, ji_v_idx, bb);
  for (auto b : bb) {
    std::cerr << std::format("ji branch boundary: {} {}\n",
                             st.get_vertex(ii_v_idx).g_v_id(),
                             st.get_vertex(b).g_v_id());
  }

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

void compute_pre_post(pst::Tree &st, tree_meta &tm) {
  std::map<pt::idx_t, pt::idx_t> &pre = tm.pre;
  std::map<pt::idx_t, pt::idx_t> &post = tm.post;

  pt::idx_t exp_size = st.vtx_count();


  // for (pt::idx_t i = 0; i < exp_size; ++i) {
  //   pre.push_back(pc::INVALID_IDX);
  //   post.push_back(pc::INVALID_IDX);
  // }

  for (pt::idx_t i = 0; i < st.vtx_count(); ++i) {
    const pst::Vertex &v = st.get_vertex(i);
    pre[v.pre_order()] = i;
    post[v.post_order()] = i;
  }

  // for (pt::idx_t i {}; i < exp_size; ++i) {

  //   if (pre[i] != pc::INVALID_IDX) {
  //     std::cerr << std::format("pre: {}  found\n", i);
  //   }
  //   if (post[i] != pc::INVALID_IDX) {
  //     std::cerr << std::format("post: {}  found\n", i);
  //     break;
  //   }
  // }
}

// compute for each vertex the number of backedges starting from a descendant of
// v to an ancestor of v
std::vector<pt::idx_t> compute_be_count(pst::Tree &st) {


  std::vector<pt::idx_t> be_count(st.vtx_count(), 0);
  std::map<pt::idx_t, std::vector<pt::idx_t>> be_map; // a temp map

  // in reverse DFS
  for (pt::idx_t v_idx {st.vtx_count()} ; v_idx-- > 0 ; ) {

    const pst::Vertex &curr_v = st.get_vertex(v_idx);

    // for each child of v, check if the backedges from the child to v
    // are backedges from a descendant of v to an ancestor of v
    for (pt::idx_t c_v_idx : st.get_children(v_idx)) {
      for (pt::idx_t be_idx : be_map[c_v_idx]) {
        const pst::BackEdge &be = st.get_backedge(be_idx);

        const pst::Vertex &v = st.get_vertex(be.get_tgt());

        if (v.pre_order() < curr_v.pre_order() && v.post_order() > curr_v.post_order()) {
          // be is a backedge from a descendant of v to an ancestor of v
          be_map[v_idx].push_back(be_idx);
        }
      }

      be_map.erase(c_v_idx);
    }

    be_count[v_idx] = be_map[v_idx].size();

    //::idx_t be_count {0};
    std::set<std::size_t> be_idxs = st.get_obe_idxs(v_idx);

    for (auto be_idx : be_idxs) {

      if (st.get_backedge(be_idx).type() != pst::be_type_e::back_edge) {
        // filter out special types of backedges
        continue;
      }

      be_map[v_idx].push_back(be_idx);
    }
  }

  return be_count;
}

// count the number of backedges from a descendant of v to an ancestor of v
// using the difference on tree
std::vector<pt::idx_t> count_brackets(pst::Tree &st, const std::vector<pt::idx_t> &B) {

  const pt::idx_t n = st.vtx_count();

  std::vector<pt::idx_t> be_count(n, 0); // how many we’ve collected so far
  std::vector<pt::idx_t> diff(n, 0);

  // sweep over backedges
  for(pt::idx_t be_idx : B){
    const pst::BackEdge &be = st.get_backedge(be_idx);

    pt::idx_t u = be.get_src();
    pt::idx_t w = be.get_tgt();

    // strict descendant: bump at parent(u)
    if (!st.is_root(u))
      diff[st.get_parent(u)] += 1;

    // strict ancestor: subtract at w itself
    diff[w] -= 1;
  }

  // Traverse vertices in reverse DFS (so children first, then parent)
  // visits every child before its parent
  for (pt::idx_t v_idx = n; v_idx-- > 0;) {
    pt::idx_t subtotal = diff[v_idx];

    for (auto c_v_idx : st.get_children(v_idx)) {
      subtotal += be_count[c_v_idx];
    }

    be_count[v_idx] = subtotal;
  }

  // print be_count
  // for (pt::idx_t i = 0; i < be_count.size(); ++i) {
  //   std::cerr << std::format("({}, {}), ", i, be_count[i]);
  // }

  return be_count;
}


/**
 * compute the flat list for backedges from a descendant of v to an ancestor of
 * v
 *
 * @param st the tree
 * @param B the backedges
 * @param off the offset table/prefix sum
 */
std::vector<pt::idx_t>
collect_backedges_by_vertex(pst::Tree &st, const std::vector<pt::idx_t> &B,
                            const std::vector<pt::idx_t> &off,
                            std::vector<pt::idx_t> &BE) {
  const pt::idx_t n = st.vtx_count();

  // build diff array
  std::vector<pt::idx_t> diff(n, 0);
  // sweep over backedges
  for (pt::idx_t be_idx : B) {
    const pst::BackEdge &be = st.get_backedge(be_idx);

    pt::idx_t u = be.get_src();
    pt::idx_t w = be.get_tgt();

    // strict descendant: bump at parent(u)
    if (!st.is_root(u))
      diff[st.get_parent(u)] += 1;

    // strict ancestor: subtract at w itself
    diff[w] -= 1;
  }

  // 4) Allocate flat storage and a little cursor per‐vertex

  std::vector<pt::idx_t> cursor(n, 0);

  // 5) Fill in each block BE[off[v] .. off[v+1]) by
  //    walking the parent‐chain from parent(u) up to w
  for (pt::idx_t be_idx : B) {
    auto &be = st.get_backedge(be_idx);

    pt::idx_t u = be.get_src();
    pt::idx_t w = be.get_tgt();

    if (st.is_root(u))
      continue;

    pt::idx_t v = st.get_parent(u);
    while (!st.is_root(v) && v != w) {
      BE[off[v] + cursor[v]] = be_idx;
      ++cursor[v];
      v = st.get_parent(v);
    }
  }

  return BE;
}

void pre_process(pst::Tree &st, tree_meta &tm) {

  // u is in the subtree of v
  // if preorder(v) < preorder(u) and postorder(v) >= postorder(u)
  // or preorder(v) < preorder(u) and preorder(u) < postorder(v) ?? confirm?

  // populate B
  std::vector<pt::idx_t> &B = tm.B;
  for (pt::idx_t be_idx {}; be_idx < st.back_edge_count(); ++be_idx) {
    const pst::BackEdge & be = st.get_backedge(be_idx);

    if (be.type() != pst::be_type_e::back_edge) {
      continue;
    }

    B.push_back(be_idx);
  }

  // compute be count
  std::vector<pt::idx_t> count = count_brackets(st, B);

  std::vector<pt::idx_t> &off = tm.off;

  // init off
  for (pt::idx_t i = 0; i < st.vtx_count()+1; ++i) {
    off.push_back(0);
  }
  //std::vector<pt::idx_t> off(st.vtx_count() + 1, 0);
  for (pt::idx_t i = 0; i < st.vtx_count(); ++i) {
    off[i+1] = off[i] + count[i];
  }

  std::vector<pt::idx_t> &BE = tm.BE;
  for (pt::idx_t i = 0; i < off[st.vtx_count()]; ++i) {
    BE.push_back(0);
  }
  //std::vector<pt::idx_t> BE(off[n]);

  //std::vector<pt::idx_t> BE =
  collect_backedges_by_vertex(st, B, off, BE);
}

void compute_depth(pst::Tree &st, tree_meta &tm) {
  std::vector<pt::idx_t> &depth = tm.depth; // assumed it is empty
  for (pt::idx_t v_idx {}; v_idx<  st.vtx_count(); v_idx++) {
    depth.push_back(pc::INVALID_IDX);
  }

  std::queue<pt::idx_t> q;
  pt::idx_t root_v_idx = st.get_root_idx();

  q.push(root_v_idx);
  depth[root_v_idx] = 0;

  while (!q.empty()) {
    pt::idx_t v_idx = q.front();
    q.pop();

    for (pt::idx_t c_v_idx : st.get_children(v_idx)) {
      depth[c_v_idx] = depth[v_idx] + 1;
      q.push(c_v_idx);
    }
  }

}

void find_hubbles(pst::Tree &st, const pvtr::Tree<pgt::flubble> &ft) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  std::cerr << std::format("{}\n", fn_name);

  std::vector<fl_in_tr> bar = find_in_st(st, ft);
  tree_meta tm;
  euler_tour(st, tm);
  compute_lo(st, tm);
  compute_pre_post(st, tm);
  pre_process(st, tm);
  compute_depth(st, tm);

  //tm.print();

  std::cerr << std::format("{} preprocessed\n", fn_name);

  /* a for ancestor, d for descendant */
  for (auto [si, ei] : bar) {

    mn_t mn = get_mn(st, tm, si, ei);

    std::cerr << "fl " << st.get_vertex(si).g_v_id() << " ~> " << st.get_vertex(ei).g_v_id() << "\n";

    std::cerr << std::format("m: {} n: {}\n", mn.m, mn.n);

    pt::idx_t si_v_id = st.get_vertex(si).g_v_id();
    //pt::idx_t boundary = with_top(st, si, ei);
    // for (auto b : with_ii(st, tm, mn, si, ei)) {
    //    std::cerr << std::format("boundary: {} {}\n", si_v_id, b);
    //  }

    for (auto b : with_ji(st, tm, mn, si, ei)) {
       std::cerr << std::format("boundary: {} {}\n", si_v_id, b);
    }

    //std::cerr << std::format("boundary: {} {}\n", si_v_id, boundary);
  }
}

} // namespace povu::hubbles
