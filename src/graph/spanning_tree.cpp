#include <algorithm>
#include <cstddef>
#include <format>
#include <iostream>
#include <limits>
#include <set>
#include <stack>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>


#include "./spanning_tree.hpp"


namespace spanning_tree {
const std::size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();

/*
 * Bracket
 * ----
 */
Bracket::Bracket():
  id(SIZE_T_MAX), recent_size_(0),  recent_class_(0), next_(nullptr), prev_(nullptr) {}
  Bracket::Bracket(BackEdge &b, std::size_t recent_size, std::size_t recent_class)
  : id(b.id()), recent_size_(recent_size), recent_class_(recent_class), next_(nullptr), prev_(nullptr) {
  b.set_bracket_ptr(this);
}

std::size_t Bracket::recent_class() const { return this->recent_class_; }
std::size_t Bracket::recent_size() const { return this->recent_size_; }
Bracket* Bracket::next() {return this->next_;}
Bracket* Bracket::prev() {return this->prev_;}

void Bracket::set_next(Bracket* p) {
  this->next_ = p;
}
void Bracket::set_prev(Bracket* p) { this->prev_ = p; }
void Bracket::set_recent_size(std::size_t s) { this->recent_size_ = s; }
void Bracket::set_recent_class(std::size_t c) { this->recent_class_ = c; }

/*
 * BracketList
 * ----
 */
  BracketList::BracketList(): first_(nullptr), last_(nullptr), size_(0) {}

// TODO: have one
Bracket& BracketList::top() { return *this->first_; }
Bracket* BracketList::first() { return this->first_; }
Bracket* BracketList::last() { return this->last_; }

std::size_t BracketList::size() const { return this->size_; };
void BracketList::set_size(std::size_t s) { this->size_ = s; };

void  BracketList::push(Bracket* b) {
  if (this->size() == 0) {
    this->last_ = b;
  }
  else { // > 0
    this
      ->first()
      ->set_next(b);
  }

  b->set_prev(this->first());

  /*if (this->size() > 0) {
    this->first()->set_next(b);
  }*/
  // this->first()->set_next(b);
  this->first_ = b;
  // TODO: should this inc size_?
}

void  BracketList::prepend(BracketList& b) {
  if (this->size() == 0 ) {

      if (this->last_ != nullptr && this->first_ != nullptr) {
      // TODO: throw proper exception
      std::cout << "BL prepend not null" << __LINE__ ;
      exit(1);
    }


  this->first_ = b.first();
  this->last_ = b.last();
  }
  else {

    this->last_->set_prev(b.first());
    this->last_ = b.last();

  }

  //this->set_size(this->size() + b.size());
}

void  BracketList::append(BracketList& b) {
  // assume this will always have size zero?
  // TODO: update b?

  if (this->size() == 0 ) {
    // assume first and last are nullptr

    if (this->last_ != nullptr && this->first_ != nullptr) {
      // TODO: throw proper exception
      std::cout << "BL append not null" << __LINE__ ;
      exit(1);
    }

    
    this->first_ = b.first();
    this->last_ = b.last();
  }
  else {
    
    this->last_->set_prev(b.first());
    this->last_ = b.last();
  }

  // this->last_ = b.first(); // TODO why does this work?
  //this->set_size(this->size() + b.size());
}

/*
 * Edge
 * ----
 */
Edge::Edge():  null_(true) {}
Edge::Edge(std::size_t id, std::size_t src, std::size_t tgt):
  id_(id), src(src), tgt(tgt), null_(false) {}

std::size_t Edge::id() const { return this->id_; }
std::size_t Edge::get_parent() const { return this->src; }
std::size_t Edge::get_child() const { return this->tgt; }
void Edge::set_class_idx(std::size_t c) { this->class_ = c; }
std::size_t Edge::get_class_idx() { return this->class_; }

/*
 * BackEdge
 * --------
 */
BackEdge::BackEdge(bool capping_be):
  b(nullptr),  class_(SIZE_T_MAX), capping_back_edge_(capping_be), null_(true)
{}

BackEdge::BackEdge(
  std::size_t id, std::size_t src, std::size_t tgt, bool capping_be):
  id_(id), src(src), tgt(tgt), b(nullptr), class_(SIZE_T_MAX),
  recent_class_(SIZE_T_MAX), recent_size_(SIZE_T_MAX),
  capping_back_edge_(capping_be),
  null_(false) {}

std::size_t BackEdge::id() const { return this->id_; }
std::size_t BackEdge::get_src() const { return this->src; }
std::size_t BackEdge::get_tgt() const { return this->tgt; }

std::size_t BackEdge::get_class() const { return this->class_; }
std::size_t BackEdge::get_recent_class() const { return this->recent_class_; }
std::size_t BackEdge::get_recent_size() const { return this->recent_size_; }

bool BackEdge::is_class_defined() const { return this->class_ != SIZE_T_MAX; }
bool BackEdge::is_capping_backedge() const { return this->capping_back_edge_; }

void BackEdge::set_class(std::size_t c) { this->class_ = c; }
void BackEdge::set_recent_class(std::size_t c) { this->recent_class_ = c;}
void BackEdge::set_recent_size(std::size_t s) { this->recent_size_ = s; }

Bracket* BackEdge::bracket_ptr() { return this->b; }
void BackEdge::set_bracket_ptr(Bracket* p) { this->b = p; }

/*
 * Vertex
 * ------
 */

// Constructor(s)
Vertex::Vertex():
  idx(0),
  parent_id(std::numeric_limits<size_t>::max()),
  hi_(std::numeric_limits<size_t>::max()),
  null_(true) {}

Vertex::Vertex(std::size_t id, std::size_t parent_id):
  idx (id),
  parent_id(parent_id),
  hi_(std::numeric_limits<size_t>::max()),
  null_(false){}


std::set<size_t> const &Vertex::get_ibe() const { return this->ibe; }
std::set<size_t> const &Vertex::get_obe() const { return this->obe; }
std::set<size_t> const& Vertex::get_children() const { return this->children; }
size_t const& Vertex::get_parent_idx() const { return this->parent_id; }

void Vertex::unset_null(){ this->null_ = false;}
void Vertex::add_obe(std::size_t id) { this->obe.insert(id); }
void Vertex::add_ibe(std::size_t id) { this->ibe.insert(id); }
void Vertex::add_child(std::size_t id) { this->children.insert(id); }
void Vertex::set_parent(std::size_t id) { this->parent_id = id; }
void Vertex::set_hi(std::size_t val) { this->hi_ = val; }

bool Vertex::is_root() const {
  return
    !this->null_ &&
    this->id() == 0 &&
    this->parent_id == std::numeric_limits<size_t>::max();
}


bool Vertex::is_leaf() const { return this->children.empty(); }
std::size_t  Vertex::hi() const { return this->hi_; }
std::size_t  Vertex::id() const { return this->idx; }
std::size_t  Vertex::parent() const { return this->parent_id; }

/*
 * Tree
 * ----
 */
// Constructor(s)
Tree::Tree() :
  nodes(std::vector<Vertex>{}),
  tree_edges(std::vector<Edge>{}),
  back_edges(std::vector<BackEdge>{}),
  bracket_lists(std::vector<BracketList>{}),
  equiv_class_count_(0) {}

Tree::Tree(std::size_t size) :
  nodes(std::vector<Vertex>{}),
  tree_edges(std::vector<Edge>{}),
  back_edges(std::vector<BackEdge>{}),
  bracket_lists(std::vector<BracketList>{}),
equiv_class_count_(0) {
  this->nodes.resize(size);
  this->bracket_lists.resize(size);
}


Vertex& Tree::get_root() { return this->nodes.at(0); }
std::size_t Tree::size() const { return this->nodes.size(); }
Vertex const &Tree::get_vertex(std::size_t vertex) const {
  return this->nodes.at(vertex);
}

std::size_t Tree::list_size(std::size_t vertex) {
  return this->bracket_lists.at(vertex).size();
}

std::size_t Tree::get_hi(std::size_t vertex) {
  return this->nodes.at(vertex).hi();
}

std::set<std::pair<std::size_t, std::size_t>>
Tree::get_children_w_id(std::size_t vertex) {
  std::set<std::pair<std::size_t, std::size_t>> res{};
  for (auto e_idx : this->nodes.at(vertex).get_children()) {
    res.insert(
      std::make_pair(
        this->tree_edges.at(e_idx).id(),
        this->tree_edges.at(e_idx).get_child()
      )
   );
  }
  return res;
}

std::vector<Edge> Tree::get_child_edges(std::size_t vertex) {
  std::vector<Edge> v{};

  for (auto e_idx : this->nodes.at(vertex).get_children()) {
    v.push_back(this->tree_edges.at(e_idx));
  }

  return v;
}

std::set<std::size_t> Tree::get_obe_idxs(std::size_t vertex) {
  std::set<std::size_t> v{};

  for (auto be_idx : this->nodes.at(vertex).get_obe()) {
    v.insert(be_idx);
  }

  return v;
}

std::set<std::size_t> Tree::get_ibe_idxs(std::size_t vertex) {
  // return this->nodes.at(vertex).get_ibe();
  std::set<std::size_t> v{};

  for (auto be_idx : this->nodes.at(vertex).get_ibe()) {
    v.insert(be_idx);
  }

  return v;
}

std::set<std::pair<std::size_t, std::size_t>>
Tree::get_obe_w_id(std::size_t vertex) {
   std::set<std::pair<std::size_t, std::size_t>> res{};
  for (auto e_idx : this->nodes.at(vertex).get_obe()) {
    res.insert(
      std::make_pair(
        this->back_edges.at(e_idx).id(),
        this->back_edges.at(e_idx).get_tgt())
    );
  }
  return res;
}

std::set<std::size_t> Tree::get_children(std::size_t vertex) {
  std::set<std::size_t> res{};
  for (auto e_idx : this->nodes.at(vertex).get_children()) {
    res.insert(this->tree_edges.at(e_idx).get_child());
  }
  return res;
}

std::set<std::size_t> Tree::get_ibe(std::size_t vertex) {
  std::set<std::size_t> res{};
  //if (this->nodes.at(vertex).is_root()) {}
  for (auto e_idx : this->nodes.at(vertex).get_ibe()) {
    res.insert(this->back_edges.at(e_idx).get_src());
  }
  return res;
}

std::set<std::size_t> Tree::get_obe(std::size_t vertex) {
  std::set<std::size_t> res{};
  for (auto e_idx : this->nodes.at(vertex).get_obe()) {
    res.insert(this->back_edges.at(e_idx).get_tgt());
  }
  return res;
}

bool Tree::is_root(std::size_t vertex) const {
  return vertex == 0;
}
bool Tree::has_child(std::size_t vertex, std::size_t child_idx)  {
  return this->get_children(vertex).count(child_idx);
}

bool Tree::has_ibe(std::size_t vertex, std::size_t qry_idx)  {
  return this->get_ibe(vertex).count(qry_idx);
}

bool Tree::has_obe(std::size_t vertex, std::size_t qry_idx)  {
  return this->get_obe(vertex).count(qry_idx);
}

std::size_t Tree::get_parent(std::size_t vertex) {
  // std::cout << vertex << "\n";
    //this->tree_edges.at(vertex).get_parent();
  auto i = this->get_vertex(vertex).get_parent_idx();
  return this->tree_edges.at(i).get_parent();
}

BackEdge& Tree::get_backedge(std::size_t backedge_idx) {
  return this->back_edges.at(backedge_idx);
}

void Tree::add_vertex(Vertex &&v) {

  if (this->size() > 0 && v.is_root()) {
    std::cerr << "Issue " << v.id() << "\n"; exit(1); // TODO: throw exception
  } else if (v.is_root()) {
    this->nodes.push_back(std::move(v));
    return;
  }

  for (auto idx{this->size()}; idx <= std::max(v.id(), v.parent()); idx++) {
    this->nodes.push_back(Vertex());
  }

  this->nodes[v.id()] = std::move(v);
  this->nodes[v.parent()].add_child(v.id());
  // expect id to be size + 1?
  //this->nodes.push_back(std::move(v));
}

void Tree::add_tree_edge(std::size_t frm, std::size_t to) {
  std::size_t edge_idx = this->tree_edges.size();
  std::size_t edge_count = edge_idx + this->back_edges.size();
  this->tree_edges.push_back(std::move(Edge(edge_count, frm, to)));

  this->nodes[frm].unset_null();
  this->nodes[to].unset_null();

  this->nodes[frm].add_child(edge_idx);
  this->nodes[to].set_parent(edge_idx);
}

std::size_t Tree::add_be(std::size_t frm, std::size_t to, bool capping_be) {
  std::size_t back_edge_idx = this->back_edges.size();
  std::size_t edge_count = back_edge_idx + this->tree_edges.size();
  this->back_edges.push_back(std::move(BackEdge(edge_count, frm, to, capping_be)));

  this->nodes[frm].add_obe(back_edge_idx);
  this->nodes[to].add_ibe(back_edge_idx);

  return back_edge_idx;
}

void Tree::set_hi(std::size_t vertex, std::size_t val) {
  this->nodes.at(vertex).set_hi(val);
}
void Tree::concat_brackets(std::size_t parent_vertex, std::size_t child_vertex) {
  //BracketList& l_p = this->bracket_lists.at(parent_vertex);
  BracketList& child_bracket = this->bracket_lists.at(child_vertex);
  BracketList& parent_bracket = this->bracket_lists.at(parent_vertex);
  parent_bracket.append(child_bracket);
  child_bracket.prepend(parent_bracket);

  std::size_t ps =  parent_bracket.size();
  std::size_t cs =  child_bracket.size();

  parent_bracket.set_size(cs + ps);
  child_bracket.set_size(cs + ps);
}

void Tree::del_bracket(std::size_t vertex, std::size_t backedge_idx) {
  BackEdge& b = this->back_edges.at(backedge_idx);

  Bracket* b_ptr = b.bracket_ptr();
  if (b_ptr == nullptr) {
    std::cout << "is null" << __LINE__;
    exit(1);
  }

  Bracket* nxt_ptr = b_ptr->next();
  Bracket* prv_ptr = b_ptr->prev();

  if (nxt_ptr != nullptr) { // why can it be null?
    nxt_ptr->set_prev(prv_ptr);
  }

  if (prv_ptr != nullptr) {
    prv_ptr->set_next(nxt_ptr);
  }
  

  delete b_ptr;
  b_ptr = nullptr;

  BracketList& bl = this->bracket_lists.at(b.get_src());
  bl.set_size(bl.size() - 1);

  // TODO: is b_ptr deleted?
}

void Tree::push(std::size_t vertex, std::size_t backege_idx) {
  std::size_t s = this->bracket_lists.at(vertex).size();
  ++s;

  Bracket* b = new Bracket(this->back_edges.at(backege_idx),
                           s,
                           this->equiv_class_count_);
  this->bracket_lists.at(vertex).push(b);

  // std::size_t s = this->bracket_lists.at(vertex).size();
  this->bracket_lists.at(vertex).set_size(s);
}

Bracket& Tree::top(std::size_t vertex) {
  return this->bracket_lists.at(vertex).top();
}

std::size_t Tree::new_class() {
  return this->equiv_class_count_++;
}

Edge& Tree::get_incoming_edge(std::size_t vertex) {
  std::size_t e_idx = this->nodes.at(vertex).get_parent_idx();
  return this->tree_edges.at(e_idx);
}

void Tree::print_dot() {
  std::cout << std::format(
    "graph G {{\n"
    "\trankdir = TB;\n"
    "\tnode[shape = circle];\n"
    "\tedge [arrowhead=vee];\n"
  );

  for (std::size_t i{}; i < this->size(); i++) {
    // Vertex const& vv = this->get_vertex(i);
    for (auto c : this->get_child_edges(i)) {
      std::cout << std::format("\t{} -- {}  [label=\"{} {}\"];\n", i, c.get_child(), c.id(), c.get_class_idx());
    }


    //for (auto c : this->get_children_w_id(i)) {
      //std::cout << std::format("\t{} -- {}  [label=\"{}\"];\n", i, c.second, c.first);
    //}

    for (auto o: this->get_obe_w_id(i)) {
      std::cout << std::format("\t{} -- {} [label=\"{}\" style=dotted];\n", i, o.second, o.first );
    }
  }

  std::cout << "}" << std::endl;
}

} // namespace spanning_tree
