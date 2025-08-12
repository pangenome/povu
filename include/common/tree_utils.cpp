#include "./tree_utils.hpp"


namespace povu::tree_utils {

// compute the branching descendants of each vertex
BranchDesc br_desc(const pst::Tree &st) {

  auto is_branching = [&](pt::idx_t v_idx) -> bool {
    return (st.get_children(v_idx).size() > 1);
  };


  // key is a vertex and values are the branching descendants of the key or leaves
  std::map<pt::idx_t, std::vector<pt::idx_t>> branch_tree;
  BranchDesc branch_desc;
  std::unordered_set<pt::idx_t> explored;

  pt::idx_t root_idx = st.get_root_idx();

  std::stack<pt::idx_t> s;
  s.push(root_idx);

  // stack of braching vertices
  std::stack<pt::idx_t> br_stack;
  std::stack<pt::idx_t> edge_stack; // idx of the edge we took could be a pair with the above

  std::unordered_set<pt::idx_t> in_br_stack;
  br_stack.push(root_idx);

  auto sort_branches = [&](pt::idx_t b) {
    pt::idx_t black_e_idx = branch_desc[b].black_e_idx;

    std::vector<pt::idx_t> &temp_sorted = branch_desc[b].sorted_br;
    temp_sorted.reserve(st.get_children(b).size());

    // add the non-black edges
    for (auto [e_idx, _] : branch_desc[b].edge_data) {
      if (e_idx != black_e_idx) {
        temp_sorted.push_back(e_idx);
      }
    }

    std::sort(temp_sorted.begin(), temp_sorted.end(),
              [&](pt::idx_t a_, pt::idx_t b_) {
                return st.get_tree_edge(a_).get_child_v_idx() > st.get_tree_edge(b_).get_child_v_idx();
              });

    if (black_e_idx != pc::INVALID_IDX) {
      temp_sorted.insert(temp_sorted.begin(), black_e_idx);
    }
  };

  while(!s.empty()) {
    pt::idx_t v_idx = s.top();

    if (pv_cmp::contains(explored, v_idx)) {

      if (v_idx == root_idx) {
        sort_branches(v_idx);
      }

      s.pop();
      continue;
    }

    const pst::Vertex &v = st.get_vertex(v_idx);

    if (v_idx != root_idx && (st.get_parent_v_idx(v_idx) == root_idx || is_branching(st.get_parent_v_idx(v_idx)))) {

      pt::idx_t prt_e_idx = v.get_parent_e_idx();

      if (st.get_parent_edge(v_idx).get_color() == pgt::color_e::black) {
        branch_desc[br_stack.top()].black_e_idx = prt_e_idx;
      }

      //std::cerr << "pushing " << prt_e_idx << " onto edge stack\n";

      edge_stack.push(prt_e_idx);
    }

    if (v.is_leaf()) { // 0 children
      // if it is a leaf...
      //branch_tree[br_stack.top()].push_back(v_idx);

      //std::cerr << "top " << edge_stack.top() << " leaf " << v_idx << "\n";

      //EdgeToBranch etb{edge_stack.top(), v_idx};
      branch_desc[br_stack.top()].edge_data.insert(std::make_pair(edge_stack.top(), v_idx));

      //br_stack.pop();
      edge_stack.pop();

      explored.insert(v_idx);
      s.pop();
    }
    else if (st.get_child_count(v_idx) == 1) { // 1 child
      std::size_t c_v_idx = *st.get_children(v_idx).begin();
      s.push(c_v_idx);
      explored.insert(v_idx);
    }
    else if (st.get_children(v_idx).size() > 1) { // 2 or more, a branching path

      if (!pv_cmp::contains(in_br_stack, v_idx)) {
        // if the vertex is branching, add it to the branch stack
        br_stack.push(v_idx);
        in_br_stack.insert(v_idx);
      }

      bool is_exp = true;
      for(auto c_v_idx : st.get_children(v_idx)) {
        if (!pv_cmp::contains(explored, c_v_idx)) {
          s.push(c_v_idx);
          is_exp = false;
          break;
        }
      }

      if (is_exp) {
        pt::idx_t b = br_stack.top();
        pt::idx_t e = edge_stack.top();

        //std::cerr << "handling br " << b << " size " << branch_desc[b].edge_data.size() << "\n";
        sort_branches(b);

#ifdef DEBUG
        assert(b == v_idx);
#endif
        br_stack.pop();
        edge_stack.pop();

        branch_desc[br_stack.top()].edge_data.insert(std::make_pair(e, b));


        explored.insert(v_idx);
      }
    }
  }

  return branch_desc;
}

/**
 * compute the flat list for backedges from a descendant of v to an ancestor of
 * v
 *
 * @param st the tree
 * @param B the backedges
 * @param off the offset table/prefix sum
 */
std::vector<pt::idx_t> collect_backedges_by_vertex(const pst::Tree &st,
                                                   const std::vector<pt::idx_t> &B,
                                                   const std::vector<pt::idx_t> &off,
                                                   std::vector<pt::idx_t> &BE) {
  const pt::idx_t n = st.vtx_count();

  // build diff array
  std::vector<pt::idx_t> diff(n, 0);
  // sweep over backedges
  for (pt::idx_t be_idx : B) {
    const pst::BackEdge &be = st.get_be(be_idx);

    pt::idx_t u = be.get_src();
    pt::idx_t w = be.get_tgt();

    // strict descendant: bump at parent(u)
    if (!st.is_root(u))
      diff[st.get_parent_v_idx(u)] += 1;

    // strict ancestor: subtract at w itself
    diff[w] -= 1;
  }

  // 4) Allocate flat storage and a little cursor per‐vertex

  std::vector<pt::idx_t> cursor(n, 0);

  // 5) Fill in each block BE[off[v] .. off[v+1]) by
  //    walking the parent‐chain from parent(u) up to w
  for (pt::idx_t be_idx : B) {
    const pst::BackEdge &be = st.get_be(be_idx);

    pt::idx_t u = be.get_src();
    pt::idx_t w = be.get_tgt();

    if (st.is_root(u))
      continue;

    pt::idx_t v = st.get_parent_v_idx(u);
    while (!st.is_root(v) && v != w) {
      BE[off[v] + cursor[v]] = be_idx;
      ++cursor[v];
      v = st.get_parent_v_idx(v);
    }
  }

  return BE;
}

/**
   *
   *
   * the max depth of a vertex reached by a backedge that starts below a
   * given backedge and ends at a vertex above this vertex
   * does not include the OBE of the vertex itself
   */

void compute_LoA(const pst::Tree &st, tree_meta &tm) {
  const std::string fn_name = pv_cmp::format("[povu::hubbles::{}]", __func__);

  const std::vector<pt::idx_t> &depth = tm.depth;

  std::vector<pt::idx_t> &loa = tm.lo;
  loa.reserve(st.vtx_count());
  for (pt::idx_t i = 0; i < st.vtx_count(); ++i) {
    loa.push_back(pc::INVALID_IDX);
  }

  // 1) write a lambda whose bool(a,b) returns true when a is “lower priority”
  auto cmp = [&](pt::idx_t a, pt::idx_t b) {
    // e.g. highest priority = largest number
    return depth[a] < depth[b];
  };


  // lower depth is higher priority
  // priority is set by value with max depth
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
      const pst::BackEdge &be = st.get_be(be_idx);

      if (be.type() != pst::be_type_e::back_edge) {
        // filter out special types of backedges
        continue;
      }

      pt::idx_t tgt_v_idx = be.get_tgt();

      pq.push(tgt_v_idx);
    }
  }
}


void compute_HiD(const pst::Tree &st, tree_meta &tm) {
  const std::string fn_name = pv_cmp::format("[povu::hubbles::{}]", __func__);

  const std::vector<pt::idx_t> &depth = tm.depth;


  std::vector<pt::idx_t> &HiD = tm.HiD;
  HiD.reserve(st.vtx_count());

  for (pt::idx_t i = 0; i < st.vtx_count(); ++i) {
    HiD.push_back(pc::INVALID_IDX);
  }

  // 1) write a lambda whose bool(a,b) returns true when a is “higher priority”
  // compare greater
  auto cmp = [&](pt::idx_t a, pt::idx_t b) {
    // e.g. highest priority = lower number
    const pst::BackEdge &be_a = st.get_be(a);
    const pst::BackEdge &be_b = st.get_be(b);

    return depth[be_a.get_src()] > depth[be_b.get_src()];
  };

  // lower depth is higher priority
  // priority is set by value with max depth
  // 2) declare the pq using decltype(cmp):
  std::priority_queue<pt::idx_t, std::vector<pt::idx_t>, decltype(cmp)> pq(cmp); // pass the lambda in the vector

  for (pt::idx_t v_idx {st.vtx_count()} ; v_idx-- > 0 ; ) {
    while (!pq.empty() && st.get_be(pq.top()).get_tgt() == v_idx) {
      pq.pop();
    }

    if (!pq.empty() && !st.is_leaf(v_idx) && !st.is_root(v_idx)) {
      HiD[v_idx] = st.get_be(pq.top()).get_src();
    }

    for (pt::idx_t be_idx : st.get_obe_idxs(v_idx)){
      pq.push(be_idx);
    }
  }
}

void compute_bracket_vals(const pst::Tree &st, tree_meta &tm) {
  const std::string fn_name = pv_cmp::format("[povu::hubbles::{}]", __func__);

  compute_LoA(st, tm);
  compute_HiD(st, tm);
}


void euler_tour(const pst::Tree &st, tree_meta &tm) {

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

    if (pv_cmp::contains(explored, v_idx)) {
      s.pop();
      continue;
    }

    const pst::Vertex &v = st.get_vertex(v_idx);

    if (!v.is_root()) {
      D[v_idx] = D[st.get_parent_v_idx(v_idx)] + 1;
    }

    if (v.is_leaf()) { // 0 children
      explored.insert(v_idx);
      s.pop();
    }
    else if (st.get_child_count(v_idx) == 1) { // 1 child
      std::size_t c_v_idx = *st.get_children(v_idx).begin();
      s.push(c_v_idx);
      explored.insert(v_idx);
    }
    else if (st.get_children(v_idx).size() > 1) { // 2 or more, a branching path
      bool is_exp = true;
      for(auto c_v_idx : st.get_children(v_idx)) {
        if (!pv_cmp::contains(explored, c_v_idx)) {
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
    if (!pv_cmp::contains(seen, v_idx)) {
      seen.insert(v_idx);
      first.push_back(i);
    }
  }

  return ;
}


pt::idx_t find_lca(const tree_meta &tm, std::vector<pt::idx_t> &vtxs) {

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
  //     std::cerr << pv_cmp::format("c {} {}, ", i, D[i]);
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
std::vector<pt::idx_t> count_brackets(const pst::Tree &st, const std::vector<pt::idx_t> &B) {

  const pt::idx_t n = st.vtx_count();

  std::vector<pt::idx_t> be_count(n, 0); // how many we’ve collected so far
  std::vector<pt::idx_t> diff(n, 0);

  // sweep over backedges
  for(pt::idx_t be_idx : B){
    const pst::BackEdge &be = st.get_be(be_idx);

    pt::idx_t u = be.get_src();
    pt::idx_t w = be.get_tgt();

    // strict descendant: bump at parent(u)
    if (!st.is_root(u))
      diff[st.get_parent_v_idx(u)] += 1;

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
  //   std::cerr << pv_cmp::format("({}, {}), ", i, be_count[i]);
  // }

  return be_count;
}


void compute_depth(const pst::Tree &st, tree_meta &tm) {
  std::vector<pt::idx_t> &depth = tm.depth; // assumed it is empty
  for (pt::idx_t v_idx{}; v_idx < st.vtx_count(); v_idx++) {
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

void compute_pre_post(const pst::Tree &st, tree_meta &tm) {
  std::map<pt::idx_t, pt::idx_t> &pre = tm.pre;
  std::map<pt::idx_t, pt::idx_t> &post = tm.post;

  //pt::idx_t exp_size = st.vtx_count();

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
  //     std::cerr << pv_cmp::format("pre: {}  found\n", i);
  //   }
  //   if (post[i] != pc::INVALID_IDX) {
  //     std::cerr << pv_cmp::format("post: {}  found\n", i);
  //     break;
  //   }
  // }
}


void pre_process(const pst::Tree &st, tree_meta &tm) {

  // u is in the subtree of v
  // if preorder(v) < preorder(u) and postorder(v) >= postorder(u)
  // or preorder(v) < preorder(u) and preorder(u) < postorder(v) ?? confirm?

  // populate B
  std::vector<pt::idx_t> &B = tm.B;
  for (pt::idx_t be_idx {}; be_idx < st.back_edge_count(); ++be_idx) {
    const pst::BackEdge & be = st.get_be(be_idx);

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

tree_meta gen_tree_meta(const pst::Tree &st ) {
  const std::string fn_name = pv_cmp::format("[povu::tree_utils::{}]", __func__);

  tree_meta tm;
  euler_tour(st, tm);
  compute_depth(st, tm);
  compute_bracket_vals(st, tm);
  compute_pre_post(st, tm);
  pre_process(st, tm);

  //tm.print();

  return tm;
}
} // namespace povu::tree_utils
