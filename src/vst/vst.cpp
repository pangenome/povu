#include "./vst.hpp"
#include "../graph/graph.hpp"
#include "../graph/spanning_tree.hpp"
#include <algorithm>
#include <cstddef>
#include <format>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>
#include <map>

#include "./vst.hpp"

namespace vst {
/*
  Constants
  ---------
*/
const std::size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();

/*
  Vertex
  ------
*/

// constructor(s)
// --------------
Vertex::Vertex(std::size_t idx) :
    idx(idx),
    cycle_equiv_class_(SIZE_T_MAX),
    parent(SIZE_T_MAX),
    children(std::set<std::size_t>{})
{}

Vertex::Vertex(
  std::size_t idx, std::size_t cycle_equiv_class, std::size_t parent) :
    idx(idx),
    cycle_equiv_class_(cycle_equiv_class),
    parent(parent),
    children(std::set<std::size_t>{})
{}


// getters
// -------
std::set<std::size_t> const& Vertex::get_children() const {
  return this->children;
}

// setters
// -------
void Vertex::add_child(std::size_t child_idx) {
  this->children.insert(child_idx);
}

/*
  VST
  ---
*/

// constructor(s)
// --------------
VST::VST(spanning_tree::Tree &t) : t(std::vector<Vertex>{}) {
  this->t.reserve(t.size());

  // a map to store the VST index of the
  // largest node with a given equivalence class
  // key is equivalence class and value is the node index
  std::map<std::size_t, std::size_t> furthest;

  // add a root node
  // ---------------

  // current spanning tree vertex
  spanning_tree::Vertex curr_span_tree_v = t.get_root();

  // current variant structure tree vertex
  Vertex curr_vst_v = Vertex(curr_span_tree_v.id());
  furthest.insert({curr_span_tree_v.id(), this->root_idx()});

  this->t.push_back(curr_vst_v); // update VST

  std::size_t current_equiv_class{};
  std::size_t parent_idx{}; // parent idx in the VST

  for (std::size_t i{1}; i < t.size(); i++){

    curr_span_tree_v = t.get_vertex(i);
    current_equiv_class = t.get_incoming_edge(i).get_class_idx();
    // TODO: not hardcode zero
    parent_idx = current_equiv_class == 0
                     ? this->root_idx()
                     : furthest.at(current_equiv_class - 1);

    curr_vst_v = Vertex(curr_span_tree_v.id(), current_equiv_class, parent_idx);

    // update the map
    // --------------

    // no need to check is value is larger or is already set because
    // we are guaranteed to always add at the max value for given equiv class
    furthest.insert({current_equiv_class, this->size()});

    // update VST
    // ----------
    this->t.push_back(curr_vst_v);
    this->t.at(parent_idx).add_child(i);
  }
}

// getters
// -------
std::set<std::size_t> const& VST::get_children(std::size_t v) const {
  return this->t.at(v).get_children();
}

std::size_t VST::size() const {return this->t.size(); }
std::size_t VST::root_idx() const {return this->root_idx_; }

// IO
// --
void VST::print_dot() {
  std::cout << std::format(
    "graph G {{\n"
    "\trankdir = TB;\n"
    "\tnode[shape = circle];\n"
    "\tedge [arrowhead=vee];\n"
  );

  for (std::size_t i{}; i < this->size(); i++) {
    for (auto c : this->get_children(i)) {
      std::cout << std::format("\t{} -- {};\n", i, c);
    }
  }

  std::cout << "}" << std::endl;
}

/**
 * Compute the equivalance class of a given vertex
 */
void handle_vertex(spanning_tree::Tree& t, std::size_t v) {
  /*
   * compute v.hi
   */
  std::set<std::size_t> ibe = t.get_ibe(v);
  std::size_t hi_0 =
    ibe.empty() ? std::numeric_limits<size_t>::max() : *ibe.begin();

  // given a node v find its child with the lowest hi value
  // (closest to root)
  // its hi value is hi_1 and the dfs num of that vertex is hi_child
  // children are empty for dummy stop node

  std::set<std::size_t> children = t.get_children(v);
  std::size_t hi_1{SIZE_T_MAX};
  std::size_t hi_child{SIZE_T_MAX};
  std::size_t hi_2{std::numeric_limits<size_t>::max()};

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
   */
  // the bracketlist was created in tree constructor
  for (auto c: children) { t.concat_brackets(v, c); }

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
  std::set<std::size_t> obe = t.get_obe_idxs(v);
  for (auto b: obe) {
    t.push(v, b);
  }

  if (hi_2 < hi_0) {
    // add a capping backedge
    std::size_t be_idx = t.add_be(v, hi_2, true);
    t.push(v, be_idx);
  }

  /*
    determine equivalance class for edge v.parent() to v
   */

  // if v is not the root of the spanning tree
  if (!t.is_root(v)) {
    spanning_tree::Bracket& b = t.top(v);

    //
    if (t.list_size(v) !=  b.recent_size()) {
      b.set_recent_size(t.list_size(v));
      b.set_recent_class(t.new_class());
    }

    //std::size_t c = t.new_class();

    spanning_tree::Edge& e = t.get_incoming_edge(v);
    //std::cout << "b->rc: " << b.recent_class() << "\n";
    e.set_class_idx(b.recent_class());

    /*check for e, b equivalance*/
    if (b.recent_size() == 1) {
      b.set_recent_class(e.get_class_idx());
    }
  }
}

/**
 * Find the cycle equivalent edges
 */
void cycle_equiv(spanning_tree::Tree &t) {
  // in reverse DFS
  for (std::size_t v{t.size() - 1}; v > 0; --v) {
    // std::cout << "v: " <<  v << "\n";
    handle_vertex(t, v);
  }
}

} // namespace vst

