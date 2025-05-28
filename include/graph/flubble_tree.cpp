#include <cassert>
#include <cstddef>
#include <iostream>
#include <limits>
#include <map>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <utility>
#include <vector>
#include <algorithm>


#include "./flubble_tree.hpp"
#include "../common/types.hpp"
#include "../common/utils.hpp"
#include "spanning_tree.hpp"


namespace povu::graph::flubble_tree {

using namespace povu::graph_types;
namespace pc = povu::constants;
namespace pt = povu::types;
namespace pvtr = povu::tree;


// orientation, id, class
struct oic_t {
  pgt::or_e orientation;
  pt::id_t id;
  pt::idx_t cls;
};


inline flubble_t forwardise(std::size_t start_id, pgt::or_e start_or, std::size_t end_id, pgt::or_e end_or) {

  if (start_or == pgt::or_e::reverse && end_or == pgt::or_e::reverse) {
    return {pgt::id_or_t{end_id, pgt::or_e::forward}, pgt::id_or_t{start_id, pgt::or_e::forward}};
  }
  else {
    return {pgt::id_or_t{start_id, start_or}, pgt::id_or_t{end_id, end_or}};
  }
}


/**
  * @brief
 */
pvtr::Tree<pvst::Vertex> construct_flubble_tree(const std::vector<oic_t> &stack_,
                                                const std::vector<pt::idx_t> &next_seen) {
  std::string fn_name = std::format("[povu::algorithms::flubble_tree::{}]", __func__);

  pvtr::Tree<pvst::Vertex> ft;
  pt::id_t root_idx {};
  pvst::Vertex root_v(root_idx, pgt::id_or_t{pc::INVALID_IDX, pgt::or_e::forward},
                                             pgt::id_or_t{pc::INVALID_IDX, pgt::or_e::forward},
                                             pvst::VertexType::dummy);
  ft.add_vertex(root_v);
  ft.set_root_idx(root_idx);

  struct ci {
    pt::idx_t cl; // class
    pt::idx_t idx; // index in stack_
  };

  std::stack<ci> s;
  std::unordered_set<pt::idx_t> in_s; // classes in s

  pt::idx_t prt_v { ft.root_idx() }; // parent vertex


  for (pt::idx_t i {}; i < stack_.size(); ++i) {

    auto [or_curr, id_curr, cl_curr] = stack_[i];

    if (id_curr == pc::INVALID_IDX) {
      continue;

    }

    // find the parent vertex, applies for non-siblings
    if (in_s.contains(cl_curr)) {

      // pop until (and including) the one whose cl equals cl_curr
      while(!s.empty()){
        auto [cl, _] = s.top();
        s.pop();
        in_s.erase(cl);

        if (cl == cl_curr){
          break;
        }
      }

      //if (prt_v != ft.root_idx()) {
        // prt_v = ft.get_parent_idx(prt_v);
      //}
    }


    if ((i + 1) < next_seen[i]) {
      auto [or_nxt, id_nxt, _] = stack_[next_seen[i]];

      if (id_nxt == pc::INVALID_IDX) {
        continue;
      }

      auto [ i, j] = forwardise(id_curr, or_curr, id_nxt, or_nxt);
      //pgt::id_or_t i = {id_curr, or_curr}; // start
      //pgt::id_or_t j = {id_nxt, or_nxt}; // end

      pvst::Vertex v(id_curr, i, j, pvst::VertexType::flubble);

      pt::idx_t v_idx = ft.add_vertex(v);
      ft.add_edge(prt_v, v_idx);
      prt_v = v_idx;
    }

    s.push( {cl_curr, i} );
    in_s.insert(cl_curr);
  }

  return ft;
}

/**
 * @brief Compute the next seen index for each equivalence class
 *
 * @param stack_ The equivalence class stack
 * @return The next seen index for each equivalence class
 */
std::vector<pt::idx_t> compute_eq_class_metadata(const std::vector<oic_t> &stack_) {

  std::vector<pt::idx_t> next_seen(stack_.size(), pc::INVALID_CLS);

  // an eq class and the index of the next time it is encountered in stack_
  std::unordered_map<pt::idx_t, pt::idx_t> next_seen_map;
  next_seen_map.reserve(stack_.size());

  for (std::size_t i {stack_.size()};  i-- > 0; ) {
    auto [or_curr, id_curr, cl_curr] = stack_[i];
    pt::idx_t next_idx = next_seen_map.contains(cl_curr) ? next_seen_map[cl_curr] : i;
    //std::cerr << std::format("i: {} \t id_curr: {} \t cl_curr: {} \t next: {}\n", i, id_curr, cl_curr, next_idx);
    next_seen[i] = next_idx;
    next_seen_map[cl_curr] = i;
  }

  return next_seen;
}

/**
 * @brief Compute the equivalence class stack
 *
 * @param t The spanning tree
 * @return The equivalence class stack
 */
std::vector<oic_t> compute_eq_class_stack2(pst::Tree &t) {
  std::string fn_name = std::format("[povu::algorithms::flubble_tree::{}]", __func__);
  const pt::idx_t EXPECTED_BLACK_EDGE_COUNT = t.vtx_count() / 2;

  pt::idx_t root_idx = t.get_root_idx();

  std::vector<oic_t> stack;
  stack.reserve(EXPECTED_BLACK_EDGE_COUNT);

  for (pt::idx_t v_idx{t.vtx_count()}; v_idx-- > 0; ) {

    if (v_idx == root_idx) { continue; }

    const pst::Edge &e = t.get_parent_edge(v_idx);
    // std::cerr << "i: " << v_idx << "\tid: " << t.get_vertex(v_idx).g_v_id() << "\tclass: " << e.get_class() << "\tclr: " << e.get_color() << "\n";

    if (e.get_color() == color::black) {
      or_e o = t.get_vertex(v_idx).type() == pgt::v_type_e::r ? pgt::or_e::forward : pgt::or_e::reverse;
      stack.push_back({o, t.get_vertex(v_idx).g_v_id(), e.get_class()});
    }
    else {
      stack.push_back({pgt::or_e::forward, pc::INVALID_IDX, e.get_class()});
    }
  }

  // print stack
  std::cerr << fn_name << " stack size: " << stack.size() << "\n";

  for (pt::idx_t i = 0; i < stack.size(); ++i) {
    const auto &oic = stack[i];
    std::cerr << "i: " << i
              << " id " << oic.id
              << " orientation " << oic.orientation
              << " class " << oic.cls << "\n";
  }

  stack.shrink_to_fit();
  std::reverse(stack.begin(), stack.end());

  return stack;
}

// compute the branching descendants of each vertex
std::map<pt::idx_t, std::vector<pt::idx_t>> br_desc(pst::Tree &t) {
  std::map<pt::idx_t, std::vector<pt::idx_t>> branch_tree;

  auto is_branching = [&](pt::idx_t v_idx) -> bool {
    return (t.get_children(v_idx).size() > 1);
  };

  std::unordered_set<pt::idx_t> explored;

  std::stack<pt::idx_t> s;
  s.push(t.get_root_idx());

  std::stack<pt::idx_t> br_stack;
  std::unordered_set<pt::idx_t> in_br_stack;
  br_stack.push(t.get_root_idx());

  while(!s.empty()) {
    pt::idx_t v_idx = s.top();

    if (explored.contains(v_idx)) {
      s.pop();
      continue;
    }

    const pst::Vertex &v = t.get_vertex(v_idx);

    if (v.is_leaf()) { // 0 children
      branch_tree[br_stack.top()].push_back(v_idx);
      explored.insert(v_idx);
      s.pop();
    }
    else if (t.get_child_edges(v_idx).size() == 1) { // 1 child
      std::size_t c_v_idx = *t.get_children(v_idx).begin();
      s.push(c_v_idx);
      explored.insert(v_idx);
    }
    else if (t.get_children(v_idx).size() > 1) { // 2 or more, a branching path

      if (!in_br_stack.contains(v_idx)) {
        // if the vertex is branching, add it to the branch stack
        br_stack.push(v_idx);
        in_br_stack.insert(v_idx);
      }

      bool is_exp = true;
      for(auto c_v_idx : t.get_children(v_idx)) {
        if (!explored.contains(c_v_idx)) {
          s.push(c_v_idx);
          is_exp = false;
          break;
        }
      }

      if (is_exp) {
        pt::idx_t b = br_stack.top();

#ifdef DEBUG
        assert(b == v_idx);
#endif
        br_stack.pop();

        branch_tree[br_stack.top()].push_back(b);

        explored.insert(v_idx);
      }
    }
  }

  return branch_tree;
}

std::vector<oic_t> compute_eq_class_stack(pst::Tree &t) {
  std::string fn_name = std::format("[povu::algorithms::flubble_tree::{}]", __func__);
  const pt::idx_t EXPECTED_BLACK_EDGE_COUNT = t.vtx_count() / 2;

  std::map<pt::idx_t, std::vector<pt::idx_t>> branch_tree = br_desc(t);

  auto is_branching = [&](pt::idx_t v_idx) -> bool {
    return (t.get_children(v_idx).size() > 1);
  };

  pt::idx_t root_idx = t.get_root_idx();

  std::vector<oic_t> stack;
  stack.reserve(EXPECTED_BLACK_EDGE_COUNT);

  // v idx to stack
  std::map<pt::idx_t, std::vector<oic_t>> stack2; // br_vtx_map

  // edge idx to mini stack
  std::map<pt::idx_t, std::vector<oic_t>> cache; // edge_map

  // e_idx to next branching vertex
  //std::map<pt::idx_t, pt::idx_t> e_idx_to_br;
  pt::idx_t last_br {pc::INVALID_IDX};

  std::vector<oic_t> mini_stack;

  for (pt::idx_t v_idx{t.vtx_count()}; v_idx-- > 0; ) {

    if (t.is_leaf(v_idx) || is_branching(v_idx)) {
      last_br = v_idx;
    }

    if (v_idx == root_idx) {

      std::vector<pt::idx_t> ce = t.get_child_edge_idxs(v_idx);

      pt::idx_t add_last{pc::INVALID_IDX};
      for (auto e_idx : ce) {
        if (t.get_tree_edge(e_idx).get_color() == color::black) {
          // by definition of the bidirected graph this happens once
          add_last = e_idx;
          // remove add_last from ce
          ce.erase(std::remove(ce.begin(), ce.end(), add_last), ce.end());
        }
      }

      std::sort(ce.begin(), ce.end(), [&](pt::idx_t a, pt::idx_t b) {
        return t.get_tree_edge(a).get_child_v_idx() <
               t.get_tree_edge(b).get_child_v_idx();
      });

      if (add_last != pc::INVALID_IDX) {
        ce.push_back(add_last);
      }



      // do in reverse
      for (auto it = ce.rbegin(); it != ce.rend(); ++it) {
        pt::idx_t ch_e_idx = *it;

        const auto &stackette = cache[ch_e_idx];

        std::cerr << "c  " << t.get_tree_edge(ch_e_idx).id() << " size ";
        std::cerr << stackette.size() << " \n";

        for (const auto &oic : stackette) {
          stack2[v_idx].push_back(oic);
        }
        // cache.erase(ch_e_idx);
      }

      continue; // TODO: is this continue needed?
    }

    const pst::Edge &e = t.get_parent_edge(v_idx);

    if (e.get_color() == color::black) {
      or_e o = t.get_vertex(v_idx).type() == pgt::v_type_e::r ? pgt::or_e::forward : pgt::or_e::reverse;
      mini_stack.push_back({o, t.get_vertex(v_idx).g_v_id(), e.get_class()});
    }
    else {
      mini_stack.push_back({pgt::or_e::forward, pc::INVALID_IDX, e.get_class()});
    }


    // if the parent is braching
    if ((!t.is_root(v_idx) && is_branching(t.get_parent(v_idx))) || (t.is_root(t.get_parent(v_idx))) ) {
      const pst::Edge &e = t.get_parent_edge(v_idx);
      pt::idx_t e_idx = t.get_vertex(v_idx).get_parent_e_idx();

      if (!t.is_leaf(last_br)) {
        cache[e_idx] = stack2[last_br];
        std::cerr << "e id " << e.id() << " lb " << last_br << "\n";
        stack2.erase(last_br);
      }

      for (const auto &oic : mini_stack) {
        cache[e_idx].push_back(oic);
      }

      mini_stack.clear();
    }

    // merge the needed edges
    // if the current vertex is branching
    if (is_branching(v_idx)) {

      std::vector<pt::idx_t> ce = t.get_child_edge_idxs(v_idx);

      pt::idx_t add_last { pc::INVALID_IDX };
      for (auto e_idx : ce) {
        if (t.get_tree_edge(e_idx).get_color() == color::black) {
          // by definition of the bidirected graph this happens once
          add_last = e_idx;
          // remove add_last from ce
          ce.erase(std::remove(ce.begin(), ce.end(), add_last), ce.end());
        }
      }

      std::sort(ce.begin(), ce.end(), [&](pt::idx_t a, pt::idx_t b) {
        return t.get_tree_edge(a).get_child_v_idx() < t.get_tree_edge(b).get_child_v_idx();
      });

      if (add_last != pc::INVALID_IDX) {
        ce.push_back(add_last);
      }

      // do in reverse
      for (auto it = ce.rbegin(); it != ce.rend(); ++it) {
        pt::idx_t ch_e_idx = *it;
        const auto &stackette = cache[ch_e_idx];

        for (const auto &oic : stackette) {
          stack2[v_idx].push_back(oic);
        }
        cache.erase(ch_e_idx);
      }
      //cache.clear();
    }

  }

  stack = stack2[root_idx];


  stack.shrink_to_fit();
  std::reverse(stack.begin(), stack.end());

  // print stack
  std::cerr << fn_name << " stack size: " << stack.size() << "\n";
  for (pt::idx_t i = 0; i < stack.size(); ++i) {
    const auto &oic = stack[i];
    std::cerr << "i: " << i
              << " id " << oic.id
              << " orientation " << oic.orientation
              << " class " << oic.cls << "\n";
  }

  return stack;
}

pvtr::Tree<pvst::Vertex> st_to_ft(pst::Tree &t) {
  std::string fn_name = std::format("[povu::algorithms::{}]", __func__);

  std::vector<oic_t> s { compute_eq_class_stack(t) };
  std::vector<pt::idx_t> next_seen = compute_eq_class_metadata(s);
  //pvtr::Tree<pvtr::Vertex<pvst::Vertex>> ft construct_flubble_tree(s, next_seen);
  pvtr::Tree<pvst::Vertex> ft = construct_flubble_tree(s, next_seen);

  return ft;
}
}
