#include "./vst.hpp"
#include "../graph/graph.hpp"
#include "../graph/spanning_tree.hpp"
#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <format>
#include <iostream>
#include <limits>
#include <utility>
#include <vector>
#include <map>
#include <stack>

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

  for (std::size_t i{1}; i < t.size(); i++){

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

} // namespace vst

namespace tree {
  // TODO: declare all of these in a header file
  const std::size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();

  // Vertex
  // ======
  
  // constructor(s)
  
  Vertex::Vertex() :
    id(SIZE_T_MAX), parent(SIZE_T_MAX), children(std::set<std::size_t>{}), is_valid_(false) {}
  Vertex::Vertex(std::size_t id) :
    id(id), parent(SIZE_T_MAX), children(std::set<std::size_t>{}), is_valid_(true) {}
  Vertex::Vertex(std::size_t id, std::size_t parent_id) :
    id(id), parent(parent_id), children(std::set<std::size_t>{}), is_valid_(true) {}

  // member function(s)
  // ------------------

  // getters
  bool Vertex::is_valid() const { return this->is_valid_; }
  std::set<std::size_t> const& Vertex::get_children() const { return this->children; }

  // setters
  // -------
  
  void Vertex::add_child(std::size_t child_id) {
    this->children.insert(child_id);
  }
  
  // Tree
  // ====

  // constructor(s)
  
  Tree::Tree() : vertices(std::vector<tree::Vertex>{}) {}

  Tree::Tree(std::size_t n) : vertices(std::vector<Vertex>(n, Vertex())) {}

  // member function(s)
  // ------------------

  // setters
  bool Tree::add_vertex(std::size_t parent_id, std::size_t id) {
    // TODO: there's a logical error in the caller if the vertex is already in the tree
    //       should we throw an exception here?
    if (this->vertices[id].is_valid()) { return false;  }
    this->vertices[id] = Vertex(id, parent_id);
    this->vertices[parent_id].add_child(id);
    return true;
  }

  std::set<std::size_t> const& Tree::get_children(std::size_t v) const {
    return this->vertices.at(v).get_children();
  }

  void Tree::print_dot() {
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

  // TODO: should this be the number of valid vertices?
  std::size_t Tree::size() const { return this->vertices.size(); }
  
} // namespace tree


namespace pst {

// constructor(s)
// --------------
PST::PST(spanning_tree::Tree &t) : t(std::vector<vst::Vertex>{}) {
  this->t.reserve(t.size()); //
  // for DFS
  std::stack<std::size_t> vertex_stack{}; // stack of vertices

  // to determine out parent child relationships in the PST
  std::stack<std::size_t> class_stack{}; // stack of classes

  // to keep track of which vertices we pushed onto the vertex stack (seen)
  std::set<std::size_t> seen{};

  // a mapping of a given class and its parent class
  std::map <std::size_t, std::size_t> class_and_parent{};

  // initialize the root class parent to itself (0)
  class_and_parent[0] = 0;

  std::size_t current_class{};
  class_stack.push(current_class);

  std::size_t current_vertex{};
  vertex_stack.push(current_vertex);

  // initialize the PST with the root vertex
  vst::Vertex root =
    vst::Vertex{current_vertex, current_class};
  this->t.push_back(root);

  vst::Vertex v = root;

  while (!vertex_stack.empty()) {
    current_vertex = vertex_stack.top();
    seen.insert(current_vertex);

    if (current_vertex > 0) {
      spanning_tree::Edge& current_edge = t.get_incoming_edge(current_vertex);
      std::size_t new_class = current_edge.get_class_idx();

      auto it = class_and_parent.find(new_class);
      if (it == class_and_parent.end()) {
        std::size_t cl_idx = this->t.size();
        class_and_parent[new_class] = cl_idx;
        //class_and_parent.find(current_class)->second;

        vst::Vertex v = vst::Vertex{current_vertex,
                                    current_class,
                                    class_and_parent.find(current_class)->second
                                   };

        this->t.push_back(v);

        if (current_class < this->t.size())   {
          this->t[current_class].add_child(cl_idx);
        } else {
          std::cout << "could not add current_vertex: " << current_vertex << "\n";
        }
      }

      current_class = new_class;
    }

    bool explored = true;

    for (auto c: t.get_children(current_vertex)){
      if (seen.find(c) == seen.end()) {
        vertex_stack.push(c);
        explored = false;
      }
    }

    if (explored) {
      vertex_stack.pop();
    }
  }
}

std::set<std::size_t> const& PST::get_children(std::size_t v) const {
  return this->t.at(v).get_children();
}

std::size_t PST::size() const {
  return this->t.size();
}

void PST::print_dot() {
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
 * PST
 * ---
 *
 * The PST is a tree of cycle equivalent regions of the spanning tree of the
 * undirected flow graph U_f with the following properties:
  * -
  *

  we use region to refer to a cycle equivalent region of the spanning tree
  therefore region is synonymous with class
  
 * The construction is as follows:
     "When a region is first entered, we set its parent to the current region 
     and then update the current region to be the region just entered.
     When a region is exited,
     the current region is set to be the exited region’s parent."
  
  * When a region is first entered:
      We know a region has been entered because: 
        - there's an OBE in the current vertex (capping or not)
    1. we set its parent to the current region
       and then update the current region to be the region just entered
       * its parent will be the class on top of the stack
         * if the stack was empty this class will be the root
    2. we will then push it onto the stack making it the current region
       the current equiv class is the class on top of the stack
   
  * When a region is exited:
      We know a region has been entered because: 
       - we see an outgoing capping back-edge
       - we see an incoming back-edge
    1 the current region is set to be the exited region’s parent.
      a) we pop the value that's on top of the stack
      b) st spanning tree  
 */
tree::Tree compute_pst(spanning_tree::Tree &st) {
  tree::Tree t = tree::Tree(st.size());
  std::stack<std::size_t> vertex_stack{}; // stack of vertices

  // go over the nodes in reverse topological order (dfs_num) which may not be
  // equal to node ids
  for (std::size_t r{st.size() - 1}; r > 0; --r) {
    // the vertex id of the vertex at the topological sort index r
    std::size_t v = st.get_sorted(r);
    
    std::set<size_t> obe = st.get_obe_idxs(v);
    std::set<size_t> ibe = st.get_ibe(v);

    if (obe.empty() && ibe.empty()) { continue; }
    
    // the equivalence class of the incoming tree edge
    std::size_t cl = st.get_parent_edge(v).get_class();

    // check if any of the obe are capping backedges
    bool has_capping_obe{false};
    for (auto e_idx: obe) {
      if (st.get_backedge(e_idx).is_capping_backedge()) {
        has_capping_obe=true;
        break;
      }
    }

    // exit the current region
    // --------------------
    
    if (!ibe.empty() || has_capping_obe) {
      vertex_stack.pop();
    }

    // enter the new region
    // --------------------

    if (!obe.empty()) { // an equivalence class/region has been entered
      if (!vertex_stack.empty()) { t.add_vertex(vertex_stack.top(), cl); }
      vertex_stack.push(cl);
    }
  }

  return t;
}
} // namespace pst


