#include "./hubbles.hpp"
#include <any>
#include <cstddef>
#include <cstdio>
#include <format>
#include <iostream>
#include <queue>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

namespace povu::hubbles {

/*
  -----------------
  Types
  -----------------
*/
struct sls {
  std::vector<pgt::id_or_t> ii_adj;
  std::vector<pgt::id_or_t> ji_adj;
};

struct fl_in_tr {
  pt::idx_t ii_idx; // parent idx
  // pgt::or_e i; // i orientation

  pt::idx_t ji_idx; // child idx
  // pgt::or_e j; // j orientation
};

/*
  -----------------
  Utils
  -----------------
*/

fl_in_tr fl_to_st_idxs(pst::Tree &st,
              std::map<pt::idx_t, std::pair<pt::idx_t, pt::idx_t>> g_id_to_idx,
              const pvst::Vertex &ft_v) {
  const std::string fn_name = std::format("[povu::hubbles::{}]", __func__);

  auto [start_id, start_or] = ft_v.get_start();
  auto [end_id, end_or] = ft_v.get_end();

  pt::idx_t x = (start_or == pgt::or_e::forward) ? g_id_to_idx[start_id].second
                                                 : g_id_to_idx[start_id].first;

  pt::idx_t y = (end_or == pgt::or_e::forward) ? g_id_to_idx[end_id].first
                                               : g_id_to_idx[end_id].second;

  fl_in_tr ij_pair = (st.get_vertex(x).dfs_num() < st.get_vertex(y).dfs_num())
                         ? fl_in_tr{x, y}
                         : fl_in_tr{y, x};

  return ij_pair;
}

/**
 * @brief populate g_id_to_idx_map with the indices of the left and right vertices
 */
void pop_id_idx(const pst::Tree &st, std::map<pt::idx_t, std::pair<pt::idx_t, pt::idx_t>> &g_id_to_idx_map) {

  //std::map<pt::idx_t, std::pair<pt::idx_t, pt::idx_t>> g_id_to_idx_map;

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

  //return g_id_to_idx_map;
}



bool is_trunk_vtx(pst::Tree &st, pt::idx_t si_v_idx, pt::idx_t ei_v_idx,
                  pt::idx_t query_v_idx) {
  if (st.get_vertex(si_v_idx).post_order() > st.get_vertex(query_v_idx).post_order() &&
      st.get_vertex(query_v_idx).post_order() > st.get_vertex(ei_v_idx).post_order()) {
    return true;
  }

  return false;
}

// get the strand of the slubble
pgt::id_or_t to_id_or(pst::Tree &st, std::pair<pt::idx_t, bool> a) {
  auto [v_idx, invert] = a;
  pgt::or_e o;

  switch (st.get_vertex(v_idx).type()) {
  case pst::v_type_e::l:
    o = pgt::or_e::forward;
    break;
  case pst::v_type_e::r:
    o = pgt::or_e::reverse;
    break;
  case pst::v_type_e::dummy:
    perror("Dummy vertex is slubble");
    exit(1);
    break;
  }

  if (invert) {
    switch (o) {
    case pgt::or_e::forward:
      o = pgt::or_e::reverse;
      break;
    case pgt::or_e::reverse:
      o = pgt::or_e::forward;
      break;
    default:
      perror("Invalid orientation");
      exit(1);
    }
  }

  pgt::id_or_t k{st.get_vertex(v_idx).g_v_id(), o};

  return k;
}

struct tree_meta {
  std::vector<pt::idx_t> E;
  std::vector<pt::idx_t> D;
  //std::map<pt::idx_t, std::vector<std::pair<pt::idx_t, pt::idx_t>>> branch_map;
  std::vector<pt::idx_t> first; // idx is v_idx value is the first time it is seen in E
  std::vector<pt::idx_t> lo; // LoA

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
   * does not include the OBE of the vertex itself
   */

void compute_lo(pst::Tree &st, tree_meta &tm) {
  const std::string fn_name = std::format("[povu::hubbles::{}]", __func__);

  const std::vector<pt::idx_t> &depth = tm.depth;

  std::vector<pt::idx_t> &loa = tm.lo;
  loa.reserve(st.vtx_count());
  for (pt::idx_t i = 0; i < st.vtx_count(); ++i) {
    loa.push_back(pc::INVALID_IDX);
  }

  // priority is set by value with max depth
  //std::priority_queue<pt::idx_t> pq;

  // 1) write a lambda whose bool(a,b) returns true when a is “lower priority”
  auto cmp = [&](pt::idx_t a, pt::idx_t b) {
    // e.g. highest priority = largest number
    return depth[a] < depth[b];
  };

  // 2) declare the pq using decltype(cmp):
  std::priority_queue<pt::idx_t, std::vector<pt::idx_t>, decltype(cmp)> pq(cmp); // pass the lambda in the vector

  for (pt::idx_t v_idx {st.vtx_count()} ; v_idx-- > 0 ; ) {

    while (!pq.empty() && pq.top() == v_idx) {
      pq.pop();
    }


    if (!pq.empty()) {
      loa[v_idx] = pq.top();
    }

    for (pt::idx_t be_idx : st.get_obe_idxs(v_idx)){
      const pst::BackEdge &be = st.get_backedge(be_idx);

      if (be.type() != pst::be_type_e::back_edge) {
        // filter out special types of backedges
        continue;
      }

      pt::idx_t tgt_v_idx = be.get_tgt();

      pq.push(tgt_v_idx);
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



struct mn_t {
  pt::idx_t m;
  pt::idx_t n;
};

pt::idx_t compute_m(pst::Tree &st, const tree_meta &tm, pt::idx_t ii_v_idx, pt::idx_t ji_v_idx) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};
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
    std::vector<pt::idx_t> p {src_v_idx, ji_v_idx};
    pt::idx_t l = find_lca(tm, p);

    if (height[l] >= height[ji_v_idx]) {
      std::cerr << fn_name << "hight[src] " << height[src_v_idx] << " h[ji] " << height[ji_v_idx] << "\n";
      std::cerr << fn_name << " skipping " << src_v_idx << "\n";
      continue;
    }

    if (height[l] > max_height.height) {
      // TODO: [c] remove type cast
      max_height = {(pt::idx_t)l, height[l]};
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

/**
 * check if a flubble can not contain a slubble
 * @param tm the tree meta data
 * @param ii_v_idx
 * @param ji_v_idx
 * @return true if the flubble can not contain the slubble, false otherwise
 */
bool cant_contain(const tree_meta &tm, pt::idx_t ii_v_idx, pt::idx_t ji_v_idx) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};


  const std::vector<pt::idx_t> &depth = tm.depth;
  const std::vector<pt::idx_t> &lo = tm.lo;

  std::cerr << fn_name << " " << lo[ji_v_idx] << "\n";

  if (depth[lo[ji_v_idx]] < depth[ii_v_idx]) {
    return true;
  }

  return false;
}

/*
  -----------------
  Handle S_i
  -----------------
*/

/**
get Ii trunk back edges
 */
pt::idx_t ii_trunk(pst::Tree &st, const tree_meta &tm, const mn_t &mn,
                   pt::idx_t ii_v_idx, pt::idx_t ji_v_idx) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  std::cerr << fn_name << "\n";

  auto [m, n] = mn;
  //std::vector<pt::idx_t> height = tm.depth;
  const std::vector<pt::idx_t> &height = tm.depth; // rename to depth

  std::cerr << fn_name << " m: " << m << " height[m] " << height[m] << " n "
            << n << " height[n]: " << height[n] << "\n";

  if (n == pc::INVALID_IDX || m == pc::INVALID_IDX || height[m] > height[n]) {
    return pc::INVALID_IDX; // invalid
  }

  // v_idx, height
  std::vector<std::pair<pt::idx_t, pt::idx_t>> x;

  for (pt::idx_t be_idx : st.get_ibe_idxs(ii_v_idx)) {
    const pst::BackEdge &be = st.get_backedge(be_idx);

    if (be.type() != pst::be_type_e::back_edge) {
      continue;
    }

    pt::idx_t be_src_v_idx = be.get_src();
    std::vector<pt::idx_t> p{be_src_v_idx, ji_v_idx};
    pt::idx_t l = find_lca(tm, p);
    if (height[l] > height[m]) {
      std::cerr << fn_name << " skipping " << be_src_v_idx << "\n";
      continue;
    }


    std::cerr << fn_name << " pushing back " << be_src_v_idx << "\n";
    x.push_back({l, height[l]});
  }

  //for (auto be_src_v_idx : st.get_ibe_src_v_idxs(ii_v_idx)) {

  //}

  // a is ancestor, d is descendant
  auto is_descendant =[&st](pt::idx_t a, pt::idx_t d)->bool{
    return st.get_vertex(a).post_order() > st.get_vertex(d).post_order() &&
      st.get_vertex(a).pre_order() < st.get_vertex(d).pre_order();
  };

  // sort by height in descending order
  std::sort(x.begin(), x.end(),
            [](const std::pair<pt::idx_t, pt::idx_t> &a,
               const std::pair<pt::idx_t, pt::idx_t> &b) {
              return a.second > b.second;
            });

  // loop from start to end of x and get the first one whose bracket source is
  // not from below m and is not ji already handled in x
  for (auto [k, h] : x) {
    bool invalid {false};
    for (auto be_idx : st.get_ibe_idxs(k)) {
      const pst::BackEdge &be = st.get_backedge(be_idx);
      pt::idx_t br_src_v_idx = be.get_src();

      if (is_descendant(ji_v_idx, br_src_v_idx)) {
        invalid = true;
        break;
      }

      //
      std::vector<pt::idx_t> p{ji_v_idx, br_src_v_idx};
      pt::idx_t l = find_lca(tm, p);

      if (height[l] > height[m]) {
        invalid = true;
        break;
      }
    }

    if (invalid) {
      continue;
    }

    for (auto be_idx : tm.get_brackets(k)) {
      const pst::BackEdge &be = st.get_backedge(be_idx);
      pt::idx_t br_src_v_idx = be.get_src();
      pt::idx_t br_tgt_v_idx = be.get_tgt();

      if (is_descendant(ji_v_idx, br_src_v_idx) && is_descendant(ii_v_idx, br_tgt_v_idx)) {
        invalid = true;
        break;
      }

      std::vector<pt::idx_t> p{k, br_src_v_idx};
      pt::idx_t l = find_lca(tm, p);

      if (height[l] > height[m]) {
        invalid = true;
        break;
      }
    }

    if (invalid) {
      continue;
    }

    return  k;

    //std::vector<pt::idx_t> s = { be_src_v_idx, ji_v_idx };
    //return find_lca(tm, s);
    //return be_src_v_idx;
    //return st.get_vertex(be_src_v_idx).g_v_id();
  }

  return pc::INVALID_IDX;
}

/** get Ii brach backedges */
void ii_branches(pst::Tree &st, const tree_meta &tm,
                 pt::idx_t ii_v_idx, pt::idx_t ji_v_idx,
                 std::vector<pt::idx_t> &bb) {

  // condition (i)
  if (st.get_children(ji_v_idx).size() < 2) {
    return;
  }

  const std::vector<pt::idx_t> &height = tm.D;
  const std::vector<pt::idx_t> &lo = tm.lo;


  // children who meet condition (i) and (ii)
  std::vector<pt::idx_t> x;
  for (pt::idx_t c_v_idx : st.get_children(ji_v_idx)) {
    if (st.get_vertex(c_v_idx).hi() == lo[c_v_idx] && st.get_vertex(c_v_idx).hi() == ii_v_idx) {
      std::vector<pt::idx_t> c_br = tm.get_brackets(c_v_idx);
      if (c_br.size() > 1) {
        x.push_back(c_v_idx);
      }
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
      bb.push_back(d);
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
        bb.push_back(min_depth.first);
      }

    }
  }

}

std::vector<std::pair<pt::idx_t, bool>>
with_ii(pst::Tree &st, const tree_meta &tm, const mn_t &mn, pt::idx_t ii_v_idx,
        pt::idx_t ji_v_idx) {

  std::vector<std::pair<pt::idx_t, bool>> res;
  pt::idx_t tb = ii_trunk(st, tm, mn, ii_v_idx, ji_v_idx);
  if (tb != pc::INVALID_IDX ) {
    res.push_back({tb, false});
    // std::cerr << std::format("ii trunk boundary: {} {}\n",
    //                          st.get_vertex(ii_v_idx).g_v_id(),
    //                          st.get_vertex(tb).g_v_id());
  }

  std::vector<pt::idx_t> bb; // branch boundaries
  ii_branches(st, tm, ii_v_idx, ji_v_idx, bb);
  for (auto b : bb) {
    res.push_back({b, false});
    // std::cerr << std::format("ii branch boundary: {} {}\n",
    //                          st.get_vertex(ii_v_idx).g_v_id(),
    //                          st.get_vertex(b).g_v_id());
  }

  return res;
}

/*
  -----------------
  Handle j_i
  -----------------
*/

pt::idx_t override_ji_trunk(pst::Tree &st, const tree_meta &tm, const mn_t &mn,
                            pt::idx_t ii_v_idx, pt::idx_t ji_v_idx) {

  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  //std::cerr << std::format("{}\n", fn_name);

  auto [m, n] = mn;
  const std::vector<pt::idx_t> &height = tm.depth; // rename to depth

  if (height[m] > height[n]) {
    return pc::INVALID_IDX; // invalid
  }

  std::cerr << fn_name << " m: " << m << " height[m] " << height[m] << " n " << n
            << " height[n]: " << height[n] << "\n";

  std::vector<pt::idx_t> case_ii_a_y;
  for (pt::idx_t c_v_idx: st.get_children(ji_v_idx)) {
    if (height [st.get_vertex(c_v_idx).hi()] < height [ii_v_idx]) {
      // when true the child is j_x
      std::cerr << fn_name << " skip child " << c_v_idx << "\n";
      continue;
    }

    std::cerr << fn_name << " c : " << c_v_idx << "\n";
    for (pt::idx_t be_idx : tm.get_brackets(c_v_idx)){
      const pst::BackEdge &be = st.get_backedge(be_idx);
      pt::idx_t src_v_idx = be.get_src();
      pt::idx_t tgt_v_idx = be.get_tgt();

      std::cerr << fn_name << " src: " << src_v_idx << " tgt: " << tgt_v_idx << "\n";

    }


    if (tm.get_brackets(c_v_idx).size() == 1) { // case (ii)
      std::vector<pt::idx_t> br = tm.get_brackets(c_v_idx);
      pt::idx_t be_idx = br.front();
      const pst::BackEdge &be = st.get_backedge(be_idx);
      pt::idx_t src_v_idx = be.get_src();
      pt::idx_t tgt_v_idx = be.get_tgt();



      if (height[m] < height[tgt_v_idx] && height[n] > height[tgt_v_idx]) {
        case_ii_a_y.push_back(tgt_v_idx); // case (ii) (a)
      }
    }
  }

  pt::idx_t min_v {n};
  // no branching path from source of only bracket to above j_i
  for (auto v_idx : case_ii_a_y) {
    bool valid {true};
    for (pt::idx_t be_idx: tm.get_brackets(v_idx)) {
      const pst::BackEdge &be = st.get_backedge(be_idx);
      pt::idx_t src_v_idx = be.get_src();
      pt::idx_t tgt_v_idx = be.get_tgt();

      if (height[src_v_idx] < height[ji_v_idx] || height[tgt_v_idx] > height[ii_v_idx]) {
        valid = false;
        break;
      }
    }

    if (valid && height[v_idx] < height[min_v]) {
      min_v = v_idx;
    }
  }

  if (min_v == n){
    return pc::INVALID_IDX;
  }
  else {
    return min_v;
  }

  // return pc::INVALID_IDX;

}

pt::idx_t ji_trunk(pst::Tree &st, const tree_meta &tm, const mn_t &mn,
                   pt::idx_t ii_v_idx, pt::idx_t ji_v_idx) {

  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  std::cerr << std::format("{}\n", fn_name);

  auto [m, n] = mn;
  const std::vector<pt::idx_t> &height = tm.depth;
  const std::vector<pt::idx_t> &lo = tm.lo;

  if (n == pc::INVALID_IDX || m == pc::INVALID_IDX || height[m] > height[n]) {
    return pc::INVALID_IDX; // invalid
  }

  pt::idx_t res = override_ji_trunk(st, tm, mn, ii_v_idx, ji_v_idx);
  if (res != pc::INVALID_IDX) {
    return res;
  }

  // v_idx, height
  std::vector<std::pair<pt::idx_t, pt::idx_t>> x;

  for (pt::idx_t be_idx : st.get_obe_idxs(ji_v_idx)){
    if (st.get_backedge(be_idx).type() != pst::be_type_e::back_edge) {
      continue;
    }
    pt::idx_t tgt_v_idx = st.get_backedge(be_idx).get_tgt();
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
    bool valid {true};

    std::cerr << fn_name << " tgt: " << tgt_v_idx << "\n";

    for (pt::idx_t be_idx : tm.get_brackets(tgt_v_idx)) {

      if (st.get_backedge(be_idx).type() != pst::be_type_e::back_edge) {
        continue;
      }

      pt::idx_t br_tgt_v_idx = st.get_backedge(be_idx).get_tgt();

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


    //std::cerr << "child" << c_v_idx << "\n";


    //std::cerr << "a\n";
    pt::idx_t be_idx = brackets.front(); // size should be 1
    //std::cerr << "b\n";
    const pst::BackEdge &be = st.get_backedge(be_idx);
    pt::idx_t src_v_idx = be.get_src();
    pt::idx_t tgt_v_idx = be.get_tgt();

    //std::cerr << src_v_idx << " - - -> " << tgt_v_idx << "\n";

    pt::idx_t curr_v_idx = src_v_idx;
    while (st.get_parent_v_idx(curr_v_idx) != ji_v_idx) {
      std::size_t p_v_idx = st.get_parent_v_idx(curr_v_idx);

      if (st.get_child_edges(p_v_idx).size() > 1) {
        //std::cerr << "pushing back (iv) g(c) " << p_v_idx << "\n";
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

std::vector<std::pair<pt::idx_t, bool>>
with_ji(pst::Tree &st, const tree_meta &tm,
                               const mn_t &mn, pt::idx_t ii_v_idx,
                               pt::idx_t ji_v_idx) {
  std::vector<std::pair<pt::idx_t, bool>> ends;

  pt::idx_t tb = ji_trunk(st, tm, mn, ii_v_idx, ji_v_idx);
  if (tb != pc::INVALID_IDX) {
    ends.push_back(std::make_pair(tb, true));
    // std::cerr << std::format("ji trunk boundary: {} {}\n",
    //                          st.get_vertex(ji_3v_idx).g_v_id(),
    //                          st.get_vertex(tb).g_v_id());
  }
  //std::cerr << std::format("ji trunk boundary: {} {}\n",
  //                         st.get_vertex(tb).g_v_id(),
  //                         st.get_vertex(ji_v_idx).g_v_id());

  std::vector<pt::idx_t> bb; // branch boundaries
  ji_branches(st, tm, ii_v_idx, ji_v_idx, bb);
  for (auto b : bb) {
    ends.push_back({b, false});
    // std::cerr << std::format("ji branch boundary: {} {}\n",
    //                          st.get_vertex(ii_v_idx).g_v_id(),
    //                          st.get_vertex(b).g_v_id());
  }

  return ends;
}



  void update_ft(pst::Tree &st, pvtr::Tree<pvst::Vertex> &ft,
               const tree_meta &tm,
               std::map<pt::idx_t, std::pair<pt::idx_t, pt::idx_t>> g_id_to_idx,
               std::map<pt::idx_t, sls> fl_map) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  const std::vector<pt::idx_t> &depth = tm.depth;

  pt::idx_t counter = ft.vtx_count();

  for (auto [ft_v_idx, slubbles] : fl_map) {
    const pvst::Vertex &ft_v = ft.get_vertex(ft_v_idx);

    // start and end (id and or) of the flubble
    pgt::id_or_t s = ft_v.get_start();
    pgt::id_or_t e = ft_v.get_end();

    std::cerr << fn_name << " flubble " << ft_v.as_str() << "\n";

    const std::vector<pt::idx_t> &ch = ft.get_children(ft_v_idx);


    if (!ch.empty()) {
      for (auto [st_idx_sl, or_sl] : slubbles.ii_adj) {
        pt::id_t sl_g_v_id = st.get_vertex(st_idx_sl).g_v_id();
        pgt::id_or_t k = {sl_g_v_id, or_sl};
        //pgt::flubble_t sl{s, k};
        pvst::Vertex v(counter, s, k, pvst::VertexType::slubble);
        //pvtr::Vertex<pgt::flubble_t> v = {counter, sl};
        ft.add_vertex(v);
        ft.add_edge(ft_v_idx, counter);
        counter++;
      }

      for (auto [st_idx_sl, or_sl] : slubbles.ji_adj) {
        pt::id_t sl_g_v_id = st.get_vertex(st_idx_sl).g_v_id();
        pgt::id_or_t k = {sl_g_v_id, or_sl};
        pvst::Vertex v(counter, k, e, pvst::VertexType::slubble);
        // pgt::flubble_t sl{k, e};
        // pvtr::Vertex<pgt::flubble_t> v = {counter, sl};
        ft.add_vertex(v);
        ft.add_edge(ft_v_idx, counter);
        counter++;
      }
    }
    else {
      // loop through ii_adj
        for (auto [st_idx_sl, or_sl] : slubbles.ii_adj) {
          pt::id_t sl_g_v_id = st.get_vertex(st_idx_sl).g_v_id();
          pgt::id_or_t k = {sl_g_v_id, or_sl};
          //pgt::flubble_t sl {s, k};
          //pvtr::Vertex<pgt::flubble_t> v = {counter, sl};
          pvst::Vertex v(counter, s, k, pvst::VertexType::slubble);
          ft.add_vertex(v);
          ft.add_edge(ft_v_idx, counter);

          for (std::size_t c : ch)  {
            const pvst::Vertex &c_v = ft.get_vertex(c);

            // if (!c_v.get_data().has_value()) {
            //   continue;
            // }
            if (c_v.get_type() != pvst::VertexType::flubble) {
              continue;
            }

            pt::idx_t c_idx = c_v.get_id();

            auto [c_ii, c_ji] = fl_to_st_idxs(st, g_id_to_idx, c_v);

            if (depth[st_idx_sl] > depth[c_ji]) {
              ft.del_edge(ft_v_idx, c_idx);
              ft.add_edge(c_idx, counter);
            }
          }

          counter++;
        }

        // loop through ii_adj
        for (auto [st_idx_sl, or_sl] : slubbles.ji_adj) {
          pt::id_t sl_g_v_id = st.get_vertex(st_idx_sl).g_v_id();
          pgt::id_or_t k = {sl_g_v_id, or_sl};
          //pgt::flubble_t sl{ k, e };
          //pvtr::Vertex<pgt::flubble_t> v = {counter, sl};
          pvst::Vertex v(counter, k, e, pvst::VertexType::slubble);
          ft.add_vertex(v);
          ft.add_edge(ft_v_idx, counter);

          for (std::size_t c : ch) {
            //const pvtr::Vertex<pgt::flubble_t> &c_v = ft.get_vertex(c);
            const pvst::Vertex &c_v = ft.get_vertex(c);

            // if (!c_v.get_data().has_value()) {
            //   continue;
            // }

            if (c_v.get_type() != pvst::VertexType::flubble) {
              continue;
            }

            pt::idx_t c_idx = c_v.get_id();

            auto [c_ii, c_ji] = fl_to_st_idxs(st, g_id_to_idx, c_v);

            if (depth[c_ii] > depth[st_idx_sl]) {
              ft.del_edge(ft_v_idx, c_idx);
              ft.add_edge(c_idx, counter);
            }
          }

          counter++;
        }
    }
  }

}



void find_hubbles(pst::Tree &st, pvtr::Tree<pvst::Vertex> &ft) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  //ft.print_dot();

  tree_meta tm;
  euler_tour(st, tm);
  compute_depth(st, tm);
  compute_lo(st, tm);
  compute_pre_post(st, tm);
  pre_process(st, tm);

  //tm.print();

  std::map<pt::idx_t, sls> flubble_map;

  std::map<pt::idx_t, std::pair<pt::idx_t, pt::idx_t>> g_id_to_idx;

  pop_id_idx(st, g_id_to_idx);

  for (pt::idx_t ft_v_idx{}; ft_v_idx < ft.vtx_count(); ft_v_idx++) {
    // const pvtr::Vertex<pgt::flubble_t> &ft_v = ft.get_vertex(ft_v_idx);

    const pvst::Vertex ft_v = ft.get_vertex(ft_v_idx);

    if (ft_v.get_type() != pvst::vt_e::flubble) {
      continue;
    }

    // if (!ft_v.get_data().has_value()) {
    //   continue;
    // }

    auto [si, ei] = fl_to_st_idxs(st, g_id_to_idx, ft_v);

    // find slubbles

    if (cant_contain(tm, si, ei)) {
      continue;
    }

    mn_t mn = get_mn(st, tm, si, ei);

    std::cerr << std::format("fl {} ~> {}\n", st.get_vertex(si).g_v_id(),
                             st.get_vertex(ei).g_v_id());

    std::cerr << std::format("si: {} ei: {}\n", si, ei);
    std::cerr << std::format("m: {} n: {}\n", mn.m, mn.n);

    pt::idx_t si_v_id = st.get_vertex(si).g_v_id();
    pt::idx_t ei_v_id = st.get_vertex(ei).g_v_id();

    std::vector<pgt::id_or_t> ii_adj;
    std::vector<pgt::id_or_t> ji_adj;

    for (auto b : with_ii(st, tm, mn, si, ei)) {
      pgt::id_or_t k = to_id_or(st, b);
      pgt::id_or_t k_ = {b.first, k.orientation};
      ii_adj.push_back(k_);
      std::cerr << std::format("boundary: {} {}\n", si_v_id, k.as_str());
    }

    for (auto b : with_ji(st, tm, mn, si, ei)) {
      pgt::id_or_t k = to_id_or(st, b);
      pgt::id_or_t k_ = {b.first, k.orientation};
      //ii_adj.push_back(k_);
      ji_adj.push_back(k_);
      std::cerr << std::format("boundary: {} {}\n", k.as_str(), ei_v_id);
    }

    if ((ii_adj.size() + ji_adj.size()) == 0) {
      std::cerr << std::format("{}: no slubbles found for flubble: {}\n", fn_name, ft_v.as_str());
      continue;
    }

    flubble_map[ft_v_idx] = sls{ii_adj, ji_adj};
  }

  std::cerr << fn_name << " flubble count " << ft.vtx_count() << "\n";
  
  update_ft(st, ft, tm, g_id_to_idx, flubble_map);

  std::cerr << fn_name << " flubble count " << ft.vtx_count() << "\n";
}

} // namespace povu::hubbles
