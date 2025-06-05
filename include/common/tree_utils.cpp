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
              [&](pt::idx_t a, pt::idx_t b) {
                return st.get_tree_edge(a).get_child_v_idx() > st.get_tree_edge(b).get_child_v_idx();
              });

    if (black_e_idx != pc::INVALID_IDX) {
      temp_sorted.insert(temp_sorted.begin(), black_e_idx);
    }
  };

  while(!s.empty()) {
    pt::idx_t v_idx = s.top();

    if (explored.contains(v_idx)) {

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
    else if (st.get_child_edge_idxs(v_idx).size() == 1) { // 1 child
      std::size_t c_v_idx = *st.get_children(v_idx).begin();
      s.push(c_v_idx);
      explored.insert(v_idx);
    }
    else if (st.get_children(v_idx).size() > 1) { // 2 or more, a branching path

      if (!in_br_stack.contains(v_idx)) {
        // if the vertex is branching, add it to the branch stack
        br_stack.push(v_idx);
        in_br_stack.insert(v_idx);
      }

      bool is_exp = true;
      for(auto c_v_idx : st.get_children(v_idx)) {
        if (!explored.contains(c_v_idx)) {
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


} // namespace povu::tree_utils
