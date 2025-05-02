#include <cstddef>
#include <iostream>
#include <limits>
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


inline flubble forwardise(std::size_t start_id, pgt::or_e start_or, std::size_t end_id, pgt::or_e end_or) {

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
pvtr::Tree<flubble> construct_flubble_tree(const std::vector<oic_t> &stack_, const std::vector<pt::idx_t> &next_seen) {
  std::string fn_name = std::format("[povu::algorithms::flubble_tree::{}]", __func__);

  pvtr::Tree<flubble> ft;

  struct ci {
    pt::idx_t cl; // class
    pt::idx_t idx; // index in stack_
  };

  std::stack<ci> s;
  std::unordered_set<pt::idx_t> in_s; // classes in s

  pt::idx_t prt_v { ft.root_idx() }; // parent vertex
  

  for (pt::idx_t i {}; i < stack_.size(); ++i) {

    auto [or_curr, id_curr, cl_curr] = stack_[i];

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

      if (prt_v != ft.root_idx()) {
        prt_v = ft.get_parent_idx(prt_v);
      }
    }

    if ((i + 1) < next_seen[i]) {
      auto [or_nxt, id_nxt, _] = stack_[next_seen[i]];

      flubble fl = forwardise(id_curr, or_curr, id_nxt, or_nxt);
      pvtr::Vertex<flubble> v(id_curr, fl);
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

  std::unordered_map<pt::idx_t, pt::idx_t> next_seen_map;
  next_seen_map.reserve(stack_.size());

  for (std::size_t i {stack_.size()};  i-- > 0; ) {
    auto [or_curr, id_curr, cl_curr] = stack_[i];
    pt::idx_t next_idx = next_seen_map.contains(cl_curr) ? next_seen_map[cl_curr] : i;
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
std::vector<oic_t> compute_eq_class_stack(pst::Tree &t) {
  std::string fn_name = std::format("[povu::algorithms::flubble_tree::{}]", __func__);
  const pt::idx_t EXPECTED_BLACK_EDGE_COUNT = t.vtx_count() / 2;

  std::vector<oic_t> stack;
  stack.reserve(EXPECTED_BLACK_EDGE_COUNT);

  // can this for loop condition be better?
  for (pt::idx_t v{t.vtx_count() - 1}; v > t.get_root_idx(); --v) {
    const pst::Edge &e = t.get_parent_edge(v);
    if (e.get_color() == color::black) {
      or_e o = t.get_vertex(v).type() == pgt::v_type_e::r ? pgt::or_e::forward : pgt::or_e::reverse;
      stack.push_back({o, t.get_vertex(v).g_v_id(), e.get_class()});
    }
    //else {
    //  stack.push_back({pgt::or_e::forward, pc::INVALID_ID, e.get_class()});
    //}
  }

  stack.shrink_to_fit();
  std::reverse(stack.begin(), stack.end());

  return stack;
}

pvtr::Tree<flubble> st_to_ft(pst::Tree& t) {
  std::string fn_name = std::format("[povu::algorithms::{}]", __func__);

  std::vector<oic_t> s { compute_eq_class_stack(t) };
  std::vector<pt::idx_t> next_seen = compute_eq_class_metadata(s);
  pvtr::Tree<flubble> ft = construct_flubble_tree(s, next_seen);

  return ft;
}

}
