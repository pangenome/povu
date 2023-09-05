#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <format>
#include <iostream>
#include <limits>
#include <unistd.h>
#include <utility>
#include <vector>
#include <map>
#include <stack>

//#include "./vst.hpp"
#include "../graph/tree.hpp"
//#include "../graph/u_graph.hpp"
#include "../graph/spanning_tree.hpp"

namespace algorithms {
// TODO: move to constants
/*
  Constants
  ---------
*/
const std::size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();
  

/**
 * Compute the equivalance class of a given vertex
 */
void handle_vertex(spanning_tree::Tree& t, std::size_t r) {

  std::size_t v = t.get_sorted(r);

  /*
   * compute v.hi
   * ------------
   */
  std::set<std::size_t> obe = t.get_obe(v);

  std::size_t hi_0 {SIZE_T_MAX};
  for (auto be: obe) {
    // if (t.get_vertex(be).dfs_num() < hi_0) { hi_0 = t.get_vertex(be).dfs_num(); }
    hi_0 = std::min(hi_0, t.get_vertex(be).dfs_num());
  }

  // given a node v find its child with the lowest hi value
  // (closest to root)
  // its hi value is hi_1 and the dfs num of that vertex is hi_child
  // children are empty for dummy stop node

  std::set<std::size_t> children = t.get_children(v);

  // a vector of pairs of hi values and children
  std::vector<std::pair<std::size_t, std::size_t>> hi_and_child{};
  hi_and_child.reserve(children.size());

  // for each child input the hi value and the child
  for (auto child: children) {
      hi_and_child.push_back({t.get_vertex(child).hi(), child});
  }

  // sort the vertex by hi value
  std::sort(hi_and_child.begin(), hi_and_child.end());

  // for each child vertex
  // the lowest hi value of all the children vertices
  std::size_t hi_1{SIZE_T_MAX};

  // TODO: urgent!! look into this
  // the child vertex whose hi value is hi_1
  std::size_t hi_child{SIZE_T_MAX};

  // if hi_and_child is not empty
  // assign hi_1 and hi_child from the first element of hi_and_child
  if (!hi_and_child.empty()) {
      hi_1 = hi_and_child[0].first;
      hi_child = hi_and_child[0].second;
  }

  //spanning_tree::Vertex& vv = ;
  t.get_vertex_mut(v).set_hi(std::min(hi_0, hi_1));

  //std::size_t n_hi = std::min(hi_0, hi_1);

  std::size_t hi_2{SIZE_T_MAX};

  // if hi_and_child has at least 2 elements
  // assign hi_2 from the second element of hi_and_child
  if (hi_and_child.size() > 1) {
       hi_2 = hi_and_child[1].first;
  }

  // sorts by the first element of the pair
  std::set<std::pair<std::size_t, std::size_t>> p; // TODO: rename p
  for (auto c: children) {
    p.insert(std::make_pair(t.get_hi(c), c));
  }

  if (!p.empty()) {
    hi_1 = p.begin()->first;
    hi_child = p.begin()->second;

    if (p.size() > 1) {
      hi_2 = (++p.begin())->first;
    }
  }

  std::size_t v_hi = std::min(hi_0, hi_1);
  t.set_hi(v, v_hi);

  /*
   * compute bracket list
   * --------------------
   */
  // the bracketlist was created in tree constructor
  for (auto c: children) {
    t.concat_bracket_lists(v, c);
  }

  // for each capping backedge TODO: add

  // pop incoming backedges
  // remove backedges we have reached the end of
  std::set<std::size_t> ibe_idxs = t.get_ibe_idxs(v);
  for (auto b: ibe_idxs) {
    t.del_bracket(v, b);

    // TODO: set backedge class ?? was id not enough?
    // do this in the del_bracket_method?
    spanning_tree::BackEdge& be= t.get_backedge(b);
    if (!be.is_capping_backedge() && !be.is_class_defined()) {
      be.set_class(t.new_class());
    }
  }

  // push outgoing backedges
  std::set<std::size_t> obe_i = t.get_obe_idxs(v);
  for (auto b: obe_i) {
    t.push(v, b);
  }

  if (hi_2 < hi_0) {
    // add a capping backedge
      std::size_t dest_v =  t.get_sorted(hi_2);
    std::size_t be_idx = t.add_be(v, dest_v, true);
    t.push(v, be_idx);
  }

  /*
   * determine equivalance class for edge v.parent() to v
   * ---------------------------------------------------
   */

  // if v is not the root of the spanning tree
  if (!t.is_root(v)) {
    spanning_tree::Bracket& b = t.top(v);

    if (t.list_size(v) !=  b.recent_size()) {
      b.set_recent_size(t.list_size(v));
      b.set_recent_class(t.new_class());
    }

    // when retreating out of a node the tree edge is labelled with
    // the class of the topmost bracket in the bracket stack
    spanning_tree::Edge& e = t.get_incoming_edge(v);
    e.set_class_idx(b.recent_class());
    if (b.is_capping()) {
      // e.set_class_idx(b.get_class());
      //e.set_class_idx(b.recent_class() - 1);
    } else {

    }
    //e.set_class_idx(b.recent_class());

    /*check for e, b equivalance*/
    if (b.recent_size() == 1) {
      std::size_t b_id = b.back_edge_id();
      spanning_tree::BackEdge& be = t.get_backedge_ref_given_id(b_id);
      be.set_class(e.get_class_idx());

      // b.set_recent_class(e.get_class_idx());
    }
  }
}

/**
 * Find the cycle equivalent edges
 * in reverse DFS
 */
void cycle_equiv(spanning_tree::Tree &t) {
  for (std::size_t r{t.size() - 1}; r > 0; --r) {
    // std::cout << "v: " << v << "\n";
    handle_vertex(t, r);
  }
}
} // namespace algorithms
