#include <algorithm>
#include <cstddef>
#include <limits>
#include <unistd.h>
#include <utility>
#include <vector>
#include <format>

#include "./cycle_equiv.hpp"
#include "../graph/spanning_tree.hpp"


namespace algorithms {
// TODO: move to constants
/*
  Constants
  ---------
*/
const std::size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();
using common::constants::UNDEFINED_SIZE_T;
using common::constants::INVALID_ID;

/**
 * Compute the equivalance class of a given vertex
 * Find the cycle equivalent edges
 * in reverse DFS
 */
void cycle_equiv(spanning_tree::Tree &t) {
    std::string fn_name = std::format("[povu::algorithms::{}]", __func__);
  // This for loop depends on an overflow
  for (std::size_t v{t.size() - 1}; v < UNDEFINED_SIZE_T; --v) {

    //std::cout << "v: " << v << " name: " << t.get_vertex(v).name() << std::endl;


    /*
     * compute v.hi
     * ------------
     */

    std::size_t hi_0 {UNDEFINED_SIZE_T};
    std::set<std::size_t> obe = t.get_obe(v);
    for (auto be: obe) {
      hi_0 = std::min(hi_0, t.get_vertex(be).dfs_num());
    }

    // given a node v find its child with the lowest hi value
    // (closest to root)
    // its hi value is hi_1 and the dfs num of that vertex is hi_child
    // children are empty for dummy stop node

    std::size_t hi_1 {UNDEFINED_SIZE_T};


    std::set<std::size_t> children = t.get_children(v);

    // a vector of pairs of hi values and children
    std::vector<std::pair<std::size_t, std::size_t>> hi_and_child{};
    hi_and_child.reserve(children.size());

    // for each child input the hi value and the child
    for (std::size_t child: children) {
      hi_and_child.push_back({t.get_vertex(child).hi(), child});
    }
    std::sort(hi_and_child.begin(), hi_and_child.end());
    if (!children.empty()) { hi_1 = hi_and_child[0].first; }

    t.get_vertex_mut(v).set_hi(std::min(hi_0, hi_1));

    // the child vertex whose hi value is hi_1
    std::size_t hi_child { UNDEFINED_SIZE_T };
    for (std::size_t child: children) {
      if (t.get_vertex(child).hi() == hi_1) {
        hi_child = child;
        break;
      }
    }

    // if hi_and_child has at least 2 elements
    // assign hi_2 from the second element of hi_and_child
    // works because hi_and_child is sorted
    std::size_t hi_2 { UNDEFINED_SIZE_T };
    for (std::size_t child: children) {
      if (child != hi_child) {
        hi_2 = t.get_vertex(child).hi();
        break;
      }
    }


    //std::cout << "\tcompute bracket list";
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
    for (std::size_t b: ibe_idxs) {
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
    for (std::size_t be_idx: obe_i) {
      t.push(v, be_idx);
    }

    if (hi_2 < hi_0 && !t.get_bracket_list(v).empty()) {
      // add a capping backedge
      std::size_t dest_v =  hi_2;
      std::size_t be_idx = t.add_be(v, dest_v, true);
      t.push(v, be_idx);
    }


    //std::cout << "\tdetermine equivalance class for edge v.parent() to v" << std::flush;
    /*
     * determine equivalance class for edge v.parent() to v
     * ---------------------------------------------------
     */

    // if v is not the root of the spanning tree
    if (!t.is_root(v)) {

      /*
        Add an articulating backedge
      */
      if (t.get_bracket_list(v).empty())  {
        std::size_t cl = t.new_class();
        while (true) {
          spanning_tree::Edge& e = t.get_incoming_edge(v);
          e.set_class_idx(cl);
          if (t.is_root(v-1) || !t.get_obe(v-1).empty()) { break; }
          --v;
        }
      }
      /*default behavior*/
      else {
        spanning_tree::Bracket& b = t.top(v);

        if (t.list_size(v) !=  b.recent_size()) {
          b.set_recent_size(t.list_size(v));
          b.set_recent_class(t.new_class());
        }

        // when retreating out of a node the tree edge is labelled with
        // the class of the topmost bracket in the bracket stack
        spanning_tree::Edge& e = t.get_incoming_edge(v);
        e.set_class_idx(b.recent_class());

        /*check for e, b equivalance*/
        if (b.recent_size() == 1) {
          std::size_t b_id = b.back_edge_id();
          spanning_tree::BackEdge& be = t.get_backedge_ref_given_id(b_id);
          be.set_class(e.get_class_idx());
        }
      }
    }

    //std::cout << "\tfinished loop" << std::endl;
  }
}

} // namespace algorithms
