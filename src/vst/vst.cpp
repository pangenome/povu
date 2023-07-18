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

#include "./vst.hpp"
#include "../graph/tree.hpp"
//#include "../graph/u_graph.hpp"
#include "../graph/spanning_tree.hpp"

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

Vertex::Vertex(std::size_t idx, std::size_t cycle_equiv_class) :
    idx(idx),
    cycle_equiv_class_(cycle_equiv_class),
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
 * VST
 * ---
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
  // TODO: use a reference
  // the dfs num of the root is always zero
  spanning_tree::Vertex curr_span_tree_v = t.get_root();

  // current variant structure tree vertex
  Vertex curr_vst_v = Vertex(curr_span_tree_v.dfs_num());
  furthest.insert({curr_span_tree_v.dfs_num(), this->root_idx()});

  this->t.push_back(curr_vst_v); // update VST

  std::size_t current_equiv_class{};
  std::size_t parent_idx{}; // parent idx in the VST

  for (std::size_t i{1}; i < t.size(); i++) {

    curr_span_tree_v = t.get_vertex(i);
    current_equiv_class = t.get_incoming_edge(i).get_class_idx();

    if (current_equiv_class == 0) {
      parent_idx = this->root_idx();
    }
    // more than one element
    else if (furthest.find(current_equiv_class - 1) == furthest.end()) {
      parent_idx = (--furthest.end())->second;
    }
    else {
      parent_idx = furthest.at(current_equiv_class-1);
    }

    curr_vst_v = Vertex(i, current_equiv_class, parent_idx);

    // update the map
    // --------------
    furthest.insert_or_assign(current_equiv_class, this->size());

    // update VST
    // ----------
    this->t.push_back(curr_vst_v);
    this->t.at(parent_idx).add_child(i);
  }
}

// getters
// -------
std::size_t VST::size() const { return this->t.size(); }
std::size_t VST::root_idx() const { return this->root_idx_; }
std::set<std::size_t> const& VST::get_children(std::size_t v) const {
  return this->t.at(v).get_children();
}

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

/*
 * Functions
 * ---------
 */

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
      b.set_recent_class(e.get_class_idx());
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

/*
 * PVST
 * ----
 *
 * The PVST is a tree of cycle equivalent regions of the spanning tree of the
 * undirected flow graph U_f with the following properties:

  we wish to arrange the nodes of our spanning tree such that the incoming tree
  edge, the parent edge, which could only be one edge, determines the
  equivalence class of that node. We want it such that nodes sharing the same
  equivalence class will be on the same level in the tree, will be children of
  the same parent in the tree. And then sibling branches or sibling sub-trees of
  the tree will be such that they are in tandem, the equivalence classes
  occurred in tandem or interspersed by the outer equivalence class.
  The parent node in a sub-tree shall represent the start of a bubble and the
  largest node in a sub-tree shall be the bubble end.

  cycle equivalent regions can only either be nested or disjoint.
  It therefore follows that bubbles can either only be nested or disjoint.
  Because of this nodes in our tree at a level n are all disjoint and nodes at a
  level n+1 are nested within bubbles of the nodes starting at n.

  a node at level N without children does not have any bubbles.
  A node at level N with children is a bubble start and the child with the
  largest DFS num or ID is a bubble end or topological sort ID.
  Nested the subtree of this node represents all the nested bubbles of this node
  N.

  -----------

  nodes in the outermost level of the flow graph are at the topmost level of the
  tree and nodes in the innermost bubbles are leaves in the tree

  the height of the tree gives us the deepest node in the graph

  a node in the graph can only belong to a single equivalence class <show why?>
  because bubbles can either be nested or disjoint

  a node's equivalence class is determined by the class of the incoming tree edge
  and not the outgoing tree edge. Which can easily be observed by the fact that
  a node can have multiple outgoing tree edges but only one incoming tree edge.
  If we assigned a node the class of the outgoing tree edge we risk assigning
  multiple equivalence classes to the same node leading to a contradiction of
  theorem xx.

  Traverse the spanning tree in top down topological order

  Going by theorem xx, we can only be in one of the following three states with
  regards to cycle equivalence, at any point while traversing the spanning tree:

   - a new cycle equivalent region has been entered
   - or as we traverse a path we are in the same cycle equivalent region
   - a cycle equiv region has been exited

   note that the entry of one cycle equiv region is not the exiting of another
   within the same node
   but rather we exit a cycle equiv region at node n and enter another at node
   n+1

   we the entry of one cycle equiv region isn't the same as the exit of another

  To build a PVST, it is therefore sufficient to only handle the three
  previously mentioned states.

  The construction is as follows:

  the current class is the class on top of the stack

  a stack of cycle equivalent regions
  the cycle equiv region on top of the stack is the current cycle equiv region

  at a node v

  1. Entry
  * When a region is first entered:
      We know a region has been entered because:
        - there's an IBE in the current vertex (capping or not)
    1. we set its parent to the furthest node in the top of the stack

       and then update the current region to be the region just entered
       * its parent will be the class on top of the stack
         * if the stack was empty this class will be the root

    2. we will then push it onto the stack making it the current region
       the current equiv class is the class on top of the stack

  2. Conserved
  * when a region has not changed that is the incoming tree edge class is same as the one on top of the stack
    * add the vertex as a child by looking up the parent node of the current class in a hash table
      (in the furthest of the previous class in the stack)
    * set that node to be the new furthest

  3. Exit
  * When a region is exited:
      We know a region has been exited because:
       - we see an OBE that is not a capping back-edge
    1 the current region is set to be the exited region’s parent.
      a) we pop the value that's on top of the stack

  Can a vertex have both incoming and outgoing backedges? not if we bi-edge it


  entry:
   - the stack is empty (start)
   - the parent has an incoming back-edge
 */
tree::Tree compute_pvst(spanning_tree::Tree &st) {

  tree::Tree t = tree::Tree(st.size(), true);

  // a map to store the VST index of the
  // largest (furthest) node with a given equivalence class
  // key is equivalence class and value is the node index
  std::map<std::size_t, std::size_t> furthest;

  // a map of a cycle equiv region to a parent node
  std::map<std::size_t, std::size_t> parent_map;

  std::stack<std::size_t> class_stack{}; // stack of classes

  std::size_t current_class{};

  // in forward topo sort order
  for (std::size_t r{0}; r < st.size(); r++) {
    // the vertex id of the vertex at the topological sort index r
    std::size_t v = st.get_sorted(r);

    //std::cout << "v: " << v << "\n";

    std::set<size_t> obe = st.get_obe_idxs(v);
    std::set<size_t> ibe = st.get_ibe_idxs(v);

    // there is at least one obe that is not a capping back edge
    bool has_non_capping_obe{false};
    for (std::size_t e_idx: obe) {
      if (!st.get_backedge(e_idx).is_capping_backedge()) {
        has_non_capping_obe=true;
        break;
      }
    }

    std::set<size_t> parent_ibe = std::set<size_t>{};
    if (!class_stack.empty()) {
      std::size_t prnt_idx = st.get_parent_edge(v).get_parent();
      parent_ibe = st.get_ibe_idxs(prnt_idx);
    }

    //if (!parent_ibe.empty() || class_stack.empty()) {}

    if (!parent_ibe.empty() || class_stack.empty()) { // entered new equiv class
      //std::cout << "state: 1\n" ;
      if (class_stack.empty()) { // this is the root node

        //
        // the start class is always zero

        // t.add_vertex(0, 0);

        class_stack.push(0);
        furthest.insert(std::pair<std::size_t, std::size_t>(0, 0));
        parent_map.insert(std::pair<std::size_t, std::size_t>(0, 0));
      }
      else {
        current_class = class_stack.top();
        //std::cout << "cc: " << current_class << "\n";
        std::size_t prnt_vtx = furthest.at(current_class);
        t.add_vertex(prnt_vtx, v);

        // the equivalence class of the incoming tree edge
        std::size_t cl = st.get_parent_edge(v).get_class();

        class_stack.push(cl);

        //furthest.insert(std::pair<std::size_t, std::size_t>(cl, v));
        parent_map.insert(std::pair<std::size_t, std::size_t>(cl, prnt_vtx));

        furthest.insert_or_assign(cl, v);
      }
    }
    else if (obe.empty()) { // the equiv class is maintained

      //std::cout << "state: 2\n" ;

      current_class = class_stack.top();

      //std::cout << "cc: " << current_class << "\n";

      // the equivalence class of the incoming tree edge
      //std::size_t cl = st.get_parent_edge(v).get_class();

      //current_class = class_stack.top();
      std::size_t prnt_vtx = parent_map.at(current_class);
      t.add_vertex(prnt_vtx, v);
      furthest.insert_or_assign(current_class, v);
    }
    else if (has_non_capping_obe) { // exit an equiv class

      //std::cout << "state: 3\n" ;
      current_class = class_stack.top();

      //std::cout << "cc: " << current_class << "\n";
      std::size_t prnt_vtx = parent_map.at(current_class);

      t.add_vertex(prnt_vtx, v);
      furthest.insert_or_assign(current_class, v);

      class_stack.pop();
    }
    else {
      ;
    }
  }

  return t;
}

/*
  A region has been entered because it's parent is zero or the parent edge is
  gray we could also check if it has an IBE

  build the tree as you walk down the spanning tree from 0 to n
 */
tree::Tree compute_pvst_grey(spanning_tree::Tree &st){
  tree::Tree t = tree::Tree(st.size(), true);

  std::set<std::size_t> seen{};
  std::stack<std::size_t> class_stack{}; // stack of classes

  // curr class to parent class
  std::map<std::size_t, std::size_t> parent_map{};
  std::map<std::size_t, std::size_t> furthest_map{};

  bool in{false};

  for (std::size_t r{1}; r < st.size(); r++) {

    // the vertex id of the vertex at the topological sort index r
    std::size_t v = st.get_sorted(r);

    std::cout << std::format("r:{} v:{}\n", r, v);

    spanning_tree::Edge const& parent_edge = st.get_parent_edge(v);
    core::color color_ = parent_edge.get_color();
    std::size_t class_ = parent_edge.get_class();
    std::size_t parent = parent_edge.get_parent();

    // for each node whose parent is node zero
    if (parent == 0) {
      parent_map.insert(std::make_pair(class_, 0));
      furthest_map.insert_or_assign(class_, v);
      t.add_vertex(0, v, class_);
      class_stack.push(class_);
      continue;
    }

    if (in) {

      //furthest_map.insert_or_assign(, v);

      std::size_t p = furthest_map.at(class_stack.top());

      parent_map.insert(std::make_pair(class_, p));
      t.add_vertex(p, v, class_);

      in = false;
      class_stack.push(class_);
      continue;
    }

    if (color_ == core::color::black) {
      furthest_map.insert_or_assign(class_, v);
      std::size_t p = class_stack.empty() ? parent_map.at(class_) : class_stack.top();
      t.add_vertex(p, v, class_);
      continue;
    }

    // entry and exit of a region
    if (color_ == core::color::gray) {

      std::size_t x =
        parent_map.end() == parent_map.find(class_) ? class_stack.top() : parent_map.at(class_);

      std::cout << "x: " << x << "\n";

      std::size_t tree_parent = furthest_map.at(x);

      if (seen.count(class_)) {
        while (!class_stack.empty() && class_stack.top() != class_) {
          class_stack.pop();
        }

        class_stack.pop();
        seen.erase(class_);

      } else {
        seen.insert(class_);
        //class_stack.push(class_);
        in = true;
      }

      t.add_vertex(tree_parent, v, class_);

      continue;
    }
  }
  return t;
}

} // namespace vst
