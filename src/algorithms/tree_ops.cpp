#include <algorithm>
#include <iostream>
#include <cstddef>
#include <map>
#include <set>
#include <string>
#include <sys/types.h>
#include <tuple>
#include <unistd.h>
#include <utility>
#include <vector>
#include <format>
#include <stack>
#include <queue>

#include "./cycle_equiv.hpp"

using core::constants::INVALID_ID;

namespace algorithms {
using namespace graph_types;

void print_branch_map(std::map<std::size_t, branch_vecs>& t) {
  std::cout << std::format(
    "graph G {{\n"
    "\trankdir = LR;\n"
    "\tnode[shape = circle];\n"
    "\tedge [arrowhead=vee];\n"
  );

  auto vec_to_string = [&](const std::vector<std::size_t>& vec) {
    std::string s;
    for (auto v : vec) {
      s += std::format("{} ", v);
    }
    return s;
  };

  for (auto [v_id, _] : t) {
    std::cout << std::format("\t{} [label = \"{}\"];\n", v_id, v_id) ;
  }

  for (auto [v_id, branches] : t) {

    for (auto br : branches.branches) {
      std::cout << std::format("\t{} -- {}  [label=\"{}\" color=\"{}\"];\n",
                               v_id, br.id, vec_to_string(br.data),
                               br.is_nesting ? "red" : "black");
    }

  }

  std::cout << "}" << std::endl;
}


// traverse a branch of the tree until either a dead end or a branching point
std::tuple<std::vector<std::size_t>, std::set<std::size_t>, std::size_t, bool>
traverse_branch(const spanning_tree::Edge& e,
                spanning_tree::Tree& t,
                const std::set<std::size_t>& seen,
                std::set<std::size_t>& faux_leaves,
                std::map<std::size_t, std::size_t>& faux_branches) {

  std::size_t eq_class {e.get_class()};
  std::size_t e_id {e.id()};

  std::vector<std::size_t> branch_ids; // ids of the vertices in this branch
  std::set<std::size_t> branch_classes; // classes seen in this branch

  // std::vector<std::size_t> eq_classes;
  if (e.get_color() == color::gray) {
    branch_ids.push_back(e_id);
    //eq_classes.push_back(e.get_class());
    branch_classes.insert(eq_class);
  }
  // eq_classes.push_back(eq_class);

  bool is_nesting {false}; // is this a nesting branch?
  if (seen.count(eq_class)) {
    //std::cerr << "eq class: " << eq_class << "\n";
    is_nesting = true; }

  std::size_t v_id = e.get_child();

  auto is_faux_branch = [&]() -> bool {
    if (v_id == t.get_root_idx()) { return false; }

    std::set<std::pair<std::size_t, std::size_t>> ibe = t.get_ibe_w_id(v_id);
    if (ibe.empty()) { return false; }

    const spanning_tree::Edge& p_e = t.get_parent_edge(v_id);

    for (auto [ibe_id, src_v_id] : ibe) {
      spanning_tree::BackEdge& be = t.get_backedge_ref_given_id(ibe_id);

      if (be.is_capping_backedge()) { continue; }

      for (auto c_e : t.get_child_edges(src_v_id)) {
        if (c_e.get_class() == p_e.get_class()) {
          faux_branches.insert({v_id, c_e.id()});
          faux_leaves.insert(src_v_id);
          return true;
        }
      }
    }

    return false;
  };

  while (true) {
    std::vector<spanning_tree::Edge> child_edges = t.get_child_edges(v_id);

    std::vector<std::size_t> be_ids;
    std::set<std::size_t> be_classes;
    for (std::size_t be_idx : t.get_obe_idxs(v_id)) {
      spanning_tree::BackEdge& be = t.get_backedge(be_idx);
      if (be.get_color() == color::gray && be.get_class() != core::constants::UNDEFINED_SIZE_T) {
        be_ids.push_back(be.id());
        be_classes.insert(be.get_class());
      }
    }

    if (faux_leaves.count(v_id) || is_faux_branch() || child_edges.size() != 1) {
      // append be_classes to eq_classes
      branch_ids.insert(branch_ids.end(), be_ids.begin(), be_ids.end());
      branch_classes.insert(be_classes.begin(), be_classes.end());

      if (child_edges.size() == 0 || faux_leaves.count(v_id)) {
        v_id = core::constants::INVALID_ID;
      }

      break;
    }

    eq_class = child_edges.front().get_class();
    e_id = child_edges.front().id();

    if (child_edges.front().get_color() == color::gray) {
      //eq_classes.push_back(eq_class);
      branch_ids.push_back(e_id);
      branch_classes.insert(eq_class);
      //if (seen.count(eq_class)) { is_nesting = true; }
    }

    branch_ids.insert(branch_ids.end(), be_ids.begin(), be_ids.end());
    branch_classes.insert(be_classes.begin(), be_classes.end());

    if (seen.count(eq_class)) {
      //std::cerr << "eq class: " << eq_class << "\n";
      is_nesting = true;
    }

    v_id = child_edges.front().get_child();
  }

  return std::make_tuple(branch_ids, branch_classes, v_id, is_nesting);
}


std::vector<std::size_t> merge_banches(std::map<std::size_t, branch_vecs>& m) {
  std::string fn_name = std::format("[povu::main::{}]", __func__);

  std::size_t root_idx {};
  std::stack<std::size_t> s;
  s.push(root_idx);

  auto is_leaf([&](std::size_t c_id) {
    //const std::vector<branch>& b = m.at(c_id).branches;
    const std::vector<branch>& branches = m.at(c_id).branches;
    return std::all_of(branches.begin(), branches.end(), [&](const branch& b) { return b.id == INVALID_ID; });
  });

  while (!s.empty()) {
    std::size_t v_id = s.top();

    // std::cout << std::format("{} v_id: {}\n", fn_name, v_id);

    branch_vecs& b = m.at(v_id);

    if (is_leaf(v_id)) {
      // merge branches putting the nesting one at the bottom
      // merge the two branches based with the nesting branch to the right

      std::size_t nesting_idx {INVALID_ID};

      for (std::size_t i{} ; i < b.branches.size(); ++i) {
       branch& br = b.branches[i];
        if (!br.is_nesting) {
          b.merged.insert(b.merged.end(), br.data.begin(), br.data.end());
          br.id = INVALID_ID;
        }
        else {
          nesting_idx = i;
        }
      }


      if (nesting_idx != INVALID_ID) {
        b.merged.insert(b.merged.end(),
                        b.branches.at(nesting_idx).data.begin(),
                        b.branches.at(nesting_idx).data.end());

        b.branches.at(nesting_idx).id = INVALID_ID;
      }

      b.branches.clear();

      for (branch& br:  m.at(b.parent_id).branches) {
        if (br.id == v_id) {
          // append b.merged to br.data
          br.data.insert(br.data.end(), b.merged.begin(), b.merged.end());
          br.id = INVALID_ID;
        }
      }

      s.pop();
    }
    else {
      // push the non leaves to the stack
      for (auto br : b.branches) {
        if (br.id != INVALID_ID) {
          s.push(br.id);
        }
      }
    }
  }

  return m.at(root_idx).merged;
}

std::vector<std::size_t> compute_eq_class_stack(spanning_tree::Tree st) {
  std::set<std::size_t> seen;
  std::set<std::size_t> faux_leaves;

  // the faux vertex idx and the id of the edge that is hidden
  std::map<std::size_t, std::size_t> faux_branches;

  // a map of LCA to branch vectors
  std::map<std::size_t, branch_vecs> branch_map;

  std::queue<std::size_t> q;
  q.push(st.get_root_idx());

  while (!q.empty()) {
    std::size_t v_id = q.front();
    q.pop();

    // std::cout << std::format("v_id: {} name: {}\n", v_id, st.get_vertex(v_id).name());

    std::vector<spanning_tree::Edge> child_edges = st.get_child_edges(v_id);
    std::vector<spanning_tree::Edge> to_traverse;

    if (child_edges.size() == 0) {
      branch_map[v_id] = {};
    }
    else if (faux_branches.find(v_id) != faux_branches.end()) {
      // this is a faux vertex
      //std::cout << "faux v_id: " << st.get_vertex(v_id).name() << " (" << v_id << ")" << std::endl;

      to_traverse.push_back(child_edges.front());
      //if (last_v != INVALID_ID) { q.push(last_v); }
      to_traverse.push_back(st.get_tree_edge_by_id(faux_branches[v_id]));
      //if (last_v2 != INVALID_ID) { q.push(last_v2); }
    }
    else {  // happens only for the root? or a branching point
      to_traverse = st.get_child_edges(v_id);
    }

    bool found_nesting {false};
    for (auto child_edge : to_traverse) {
      auto [e_ids, eq_cls, last_v, is_nesting] =
        traverse_branch(child_edge, st, seen, faux_leaves, faux_branches);

      if (is_nesting && !found_nesting) {
        found_nesting = true;
      }
      else if (is_nesting && found_nesting) {
        throw std::logic_error(std::format("more than 2 branches of  id {} name {} are nesting",
                                           v_id, st.get_vertex(v_id).name()) );
      }



      branch_map[v_id].branches.push_back({last_v, e_ids, is_nesting});
      seen.insert(eq_cls.begin(), eq_cls.end());

      if (!faux_leaves.count(last_v) && last_v != INVALID_ID ) {
         branch_map[last_v].parent_id = v_id;
        q.push(last_v);
      }
      else {
        branch_map[last_v] = {};
        branch_map[last_v].parent_id = v_id;
      }
    }
  }

  if (false) {
    std::cout << "branch_map\n";
    print_branch_map(branch_map);
  }

  return merge_banches(branch_map);

}


} // namespace algorithms
