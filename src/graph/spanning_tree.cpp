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
  back_edge_id_(SIZE_T_MAX),
  recent_size_(0),
  recent_class_(0) {}

Bracket::Bracket(
  std::size_t backedge_id, std::size_t recent_size, std::size_t recent_class, bool is_capping)
    : back_edge_id_(backedge_id),
      recent_size_(recent_size),
      recent_class_(recent_class),
      is_capping_(is_capping){}


bool Bracket::is_capping() const { return this->is_capping_; }  
std::size_t Bracket::back_edge_id() { return this->back_edge_id_; }
std::size_t Bracket::recent_class() const { return this->recent_class_; }
std::size_t Bracket::recent_size() const { return this->recent_size_; }
void Bracket::set_recent_size(std::size_t s) { this->recent_size_ = s; }
void Bracket::set_recent_class(std::size_t c) { this->recent_class_ = c; }

/*
 * Edge
 * ----
 */
// constructor(s)
Edge::Edge(): null_(true) {}
Edge::Edge(std::size_t id, std::size_t src, std::size_t tgt):
  id_(id), src(src), tgt(tgt), null_(false) {}

// getters
std::size_t Edge::id() const { return this->id_; }
std::size_t Edge::get_parent() const { return this->src; }
std::size_t Edge::get_child() const { return this->tgt; }
std::size_t Edge::get_class() const { return this->class_; }
std::size_t Edge::get_class_idx() { return this->class_; }

// setters
void Edge::set_class_idx(std::size_t c) { this->class_ = c; }

/*
 * BackEdge
 * --------
 */
BackEdge::BackEdge(bool capping_be):
  class_(SIZE_T_MAX), capping_back_edge_(capping_be), null_(true)
{}

BackEdge::BackEdge(std::size_t id, std::size_t src, std::size_t tgt,
                   bool capping_be)
    : id_(id), src(src), tgt(tgt), class_(SIZE_T_MAX),
      recent_class_(SIZE_T_MAX), recent_size_(SIZE_T_MAX),
      capping_back_edge_(capping_be), null_(false) {}

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
std::list<Bracket>::iterator BackEdge::bracket_it() { return this->bi; }
void BackEdge::set_bracket_it(std::list<Bracket>::iterator it) {
  this->bi = it;
}

/*
 * Vertex
 * ------
 */
// Constructor(s)
Vertex::Vertex():
  dfs_num_(SIZE_T_MAX),
  parent_id(std::numeric_limits<size_t>::max()),
  hi_(std::numeric_limits<size_t>::max()),
  null_(true) {}

Vertex::Vertex(std::size_t id, std::size_t parent_id):
  dfs_num_ (id),
  parent_id(parent_id),
  hi_(std::numeric_limits<size_t>::max()),
  null_(false){}

// getters
std::size_t  Vertex::hi() const { return this->hi_; }
std::size_t  Vertex::dfs_num() const { return this->dfs_num_; }
bool Vertex::is_leaf() const { return this->children.empty(); }
std::size_t  Vertex::parent() const { return this->parent_id; }
std::set<size_t> const &Vertex::get_ibe() const { return this->ibe; }
std::set<size_t> const &Vertex::get_obe() const { return this->obe; }
size_t const& Vertex::get_parent_idx() const { return this->parent_id; }
std::set<size_t> const& Vertex::get_children() const { return this->children; }
bool Vertex::is_root() const {
  return
    !this->null_ &&
    this->dfs_num() == 0 &&
    this->parent_id == std::numeric_limits<size_t>::max();
}

// setters    
void Vertex::unset_null(){ this->null_ = false;}
void Vertex::add_obe(std::size_t obe_id) { this->obe.insert(obe_id); }
void Vertex::add_ibe(std::size_t ibe_id) { this->ibe.insert(ibe_id); }
void Vertex::add_child(std::size_t e_id) { this->children.insert(e_id); }
void Vertex::set_parent(std::size_t n_id) { this->parent_id = n_id; }
void Vertex::set_hi(std::size_t val) { this->hi_ = val; }
void Vertex::set_dfs_num(std::size_t idx) { this->dfs_num_ = idx; }

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
  sort_(std::vector<std::size_t>{}),
  equiv_class_count_(0) {}

Tree::Tree(std::size_t size) :
  nodes(std::vector<Vertex>{}),
  tree_edges(std::vector<Edge>{}),
  back_edges(std::vector<BackEdge>{}),
  bracket_lists(std::vector<BracketList>{}),
  sort_(std::vector<std::size_t>{}),
  equiv_class_count_(0) {
  this->nodes.resize(size);
  this->bracket_lists.resize(size);
  this->sort_.resize(size);
}

void Tree::set_sort(std::size_t idx, std::size_t vertex) {
  this->sort_.at(idx) = vertex;
}
void Tree::set_dfs_num(std::size_t vertex, std::size_t dfs_num) {
  this->nodes.at(vertex).set_dfs_num(dfs_num);
}

    
Vertex& Tree::get_root() { return this->nodes.at(0); }
std::size_t Tree::size() const { return this->nodes.size(); }
Vertex const &Tree::get_vertex(std::size_t vertex) const {
  return this->nodes.at(vertex);
}
Vertex& Tree::get_vertex_mut(std::size_t vertex) {
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

Edge const& Tree::get_parent_edge(std::size_t vertex) const {
  return this->tree_edges.at(this->nodes.at(vertex).get_parent_idx());
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

BackEdge Tree::get_backedge_given_id(std::size_t backedge_id) {
  BackEdge b;
  for (auto be : this->back_edges) {
    if (be.id() == backedge_id) {
      b = be;
    }
  }
  return b;
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

/**
 * insert the elements of the child bracket list at the
 * beginning of the parent bracket list
 *
 * TODO: make constant time
 *
 * @param parent_vertex
 * @param child_vertex
*/
void
Tree::concat_bracket_lists(std::size_t parent_vertex, std::size_t child_vertex) {
  BracketList& bl_p = this->bracket_lists.at(parent_vertex);
  BracketList& bl_c = this->bracket_lists.at(child_vertex);

  bl_p.insert(bl_p.begin(), bl_c.begin(), bl_c.end() );
}

// TODO: once deleted do we care to reflect changes in the concated ones?
// WE DO!!
/**
   delete the bracket that is associated with the backedge
   given a vertex id and a backedge idx
   */
void Tree::del_bracket(std::size_t vertex, std::size_t backedge_idx) {
  BackEdge& b = this->back_edges.at(backedge_idx);

  std::list<Bracket>::iterator b_it =
    this->back_edges.at(backedge_idx).bracket_it();
    // BracketList& bl = this->bracket_lists.at(vertex);

  // auto res = this->bracket_lists.at(vertex).erase(b_it);

  for (auto br = this->bracket_lists.at(vertex).begin();
       br != this->bracket_lists.at(vertex).end();
       br++) {
    if (br->back_edge_id() == b.id()) {
      this->bracket_lists.at(vertex).erase(br);
      std::size_t s =  this->bracket_lists.at(vertex).size();
      this->bracket_lists.at(vertex).front().set_recent_size(s);
      return;
    }
  }
}

BracketList& Tree::get_bracket_list(std::size_t vertex) {
    return this->bracket_lists.at(vertex);
}

void Tree::push(std::size_t vertex, std::size_t backege_idx) {
  std::size_t s = this->bracket_lists.at(vertex).size();
  // ++s;

  BackEdge& be = this->back_edges.at(backege_idx);

  // TODO: std::move
  // back edge id, recent size, recent class
  
  Bracket br =
    Bracket(be.id(), s, this->equiv_class_count_, be.is_capping_backedge());
  this->bracket_lists.at(vertex).push_front(br);
  std::list<Bracket>::iterator br_it = this->bracket_lists.at(vertex).begin();

  // std::size_t s = this->bracket_lists.at(vertex).size();
  // this->bracket_lists.at(vertex).set_size(++s);

  // update the backedge
  be.set_bracket_it(br_it);
  // TODO: unnecssary?
  be.set_class(this->equiv_class_count_);
  be.set_recent_size(s);
}

Bracket& Tree::top(std::size_t vertex) {
  return this->bracket_lists.at(vertex).front();
}

std::size_t Tree::new_class() { return this->equiv_class_count_++; }

Edge& Tree::get_incoming_edge(std::size_t vertex) {
  std::size_t e_idx = this->nodes.at(vertex).get_parent_idx();
  return this->tree_edges.at(e_idx);
}
std::size_t Tree::get_sorted(std::size_t idx) { return  this->sort_.at(idx);}

void Tree::print_dot() {
  std::cout << std::format(
    "graph G {{\n"
    "\trankdir = TB;\n"
    "\tnode[shape = circle];\n"
    "\tedge [arrowhead=vee];\n"
  );

  for (std::size_t j{}; j < this->size(); j++){
    std::size_t i = this->get_sorted(j);
    std::cout << std::format("\t{} [label = \"{} {}_s\"];\n", i, i, j) ;
  }


  for (std::size_t j{}; j < this->size(); j++) {

    std::size_t i = this->sort_.at(j);

    // Vertex const& vv = this->get_vertex(i);
    for (auto c : this->get_child_edges(i)) {
      std::string cl = c.get_class_idx() > 10000
                           ? "\u2205"
                           : std::to_string(c.get_class_idx());

      std::cout << std::format("\t{}  -- {}  [label=\"{} {}\"];\n",
                               i, c.get_child(), c.id(), cl);
    }


    //for (auto c : this->get_children_w_id(i)) {
      //std::cout << std::format("\t{} -- {}  [label=\"{}\"];\n", i, c.second, c.first);
    //}

    for (auto o: this->get_obe_w_id(i)) {
      std::string cl = o.first > 10000 ? "\u2205" : std::to_string(o.first);
      std::string color =  this->get_backedge_given_id(o.first).is_capping_backedge() ?
        "red" :
        "black";
      
      std::cout << std::format(
          "\t{} -- {} [label=\"{}\" style=\"dotted\" color=\"{}\"];\n",
          i, o.second, cl, color
          );
    }
  }

  std::cout << "}" << std::endl;
}

} // namespace spanning_tree
