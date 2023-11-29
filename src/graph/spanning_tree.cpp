#include <algorithm>
#include <cstddef>
#include <exception>
#include <format>
#include <iostream>
#include <limits>
#include <set>
#include <stack>
#include <stdexcept>
#include <string>
#include <sys/types.h>
#include <tuple>
#include <type_traits>
#include <utility>
#include <vector>
#include <map>

#include "../core/core.hpp"
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
Edge::Edge(std::size_t id, std::size_t src, std::size_t tgt, core::color c):
  id_(id), src(src), tgt(tgt), null_(false), color_(c) {}

// getters
std::size_t Edge::id() const { return this->id_; }
core::color Edge::get_color() const { return this->color_; }
std::size_t Edge::get_parent() const { return this->src; }
std::size_t Edge::get_child() const { return this->tgt; }
std::size_t Edge::get_class() const { return this->class_; }
std::size_t Edge::get_class_idx() { return this->class_; }

// setters
void Edge::set_class_idx(std::size_t c) { this->class_ = c; }
void Edge::set_class(std::size_t c) { this->class_ = c; }

/*
 * BackEdge
 * --------
 */
BackEdge::BackEdge(bool capping_be):
  class_(SIZE_T_MAX), capping_back_edge_(capping_be), null_(true)
{}

BackEdge::BackEdge(
  std::size_t id,
  std::size_t src,
  std::size_t tgt,
  bool capping_be,
  core::color c
  )
  : id_(id), src(src), tgt(tgt), class_(SIZE_T_MAX), recent_class_(SIZE_T_MAX),
	recent_size_(SIZE_T_MAX), capping_back_edge_(capping_be), null_(false),
	color_(c) {}

std::size_t BackEdge::id() const { return this->id_; }
std::size_t BackEdge::get_src() const { return this->src; }
std::size_t BackEdge::get_tgt() const { return this->tgt; }

std::size_t BackEdge::get_class() const { return this->class_; }
std::size_t BackEdge::get_recent_class() const { return this->recent_class_; }
std::size_t BackEdge::get_recent_size() const { return this->recent_size_; }

core::color BackEdge::get_color() const { return this->color_; }

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
  sort_g(std::vector<std::size_t>{}),
  equiv_class_count_(0) {}

Tree::Tree(std::size_t size) :
  nodes(std::vector<Vertex>{}),
  tree_edges(std::vector<Edge>{}),
  back_edges(std::vector<BackEdge>{}),
  bracket_lists(std::vector<BracketList>{}),
  sort_(std::vector<std::size_t>{}),
  sort_g(std::vector<std::size_t>{}),
  equiv_class_count_(0) {
  this->nodes.resize(size);
  this->bracket_lists.resize(size);
  this->sort_.resize(size);
  this->sort_g.resize(size);
}

void Tree::set_sort(std::size_t idx, std::size_t vertex) {
  this->sort_.at(idx) = vertex;
}
  
void Tree::set_sort_g(std::size_t idx, std::size_t vertex) {
  this->sort_g.at(idx) = vertex;
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

BackEdge& Tree::get_backedge_ref_given_id(std::size_t backedge_id) {
  for (auto it {this->back_edges.begin()}; it != this->back_edges.end(); ++it) {
	if (it->id() == backedge_id) {
	  return *it;
	}
  }

  throw std::out_of_range(std::format("Back edge {} not found", backedge_id));
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

void Tree::add_tree_edge(std::size_t frm, std::size_t to, core::color c) {
  std::size_t edge_idx = this->tree_edges.size();
  std::size_t edge_count = edge_idx + this->back_edges.size();
  this->tree_edges.push_back(std::move(Edge(edge_count, frm, to, c)));

  this->nodes[frm].unset_null();
  this->nodes[to].unset_null();

  this->nodes[frm].add_child(edge_idx);
  this->nodes[to].set_parent(edge_idx);
}

std::size_t Tree::add_tree_edge(std::size_t frm, std::size_t to, std::size_t weight, core::color c) {
  std::size_t edge_idx = this->tree_edges.size();
  std::size_t edge_count = edge_idx + this->back_edges.size();
  this->tree_edges.push_back(std::move(Edge(edge_count, frm, to, c)));

  this->tree_edges.at(edge_idx).set_class(weight);

  this->nodes[frm].unset_null();
  this->nodes[to].unset_null();

  this->nodes[frm].add_child(edge_idx);
  this->nodes[to].set_parent(edge_idx);

  return edge_idx;
}

std::size_t Tree::add_be(std::size_t frm, std::size_t to, bool capping_be, core::color c) {
  std::size_t back_edge_idx = this->back_edges.size();
  std::size_t edge_count = back_edge_idx + this->tree_edges.size();
  this->back_edges.push_back(std::move(BackEdge(edge_count, frm, to, capping_be, c)));

  this->nodes[frm].add_obe(back_edge_idx);
  this->nodes[to].add_ibe(back_edge_idx);

  return back_edge_idx;
}

  std::size_t Tree::add_be(std::size_t frm, std::size_t to, std::size_t weight, bool capping_be, core::color c) {
  std::size_t back_edge_idx = this->back_edges.size();
  std::size_t edge_count = back_edge_idx + this->tree_edges.size();
  this->back_edges.push_back(std::move(BackEdge(edge_count, frm, to, capping_be, c)));

  this->back_edges.at(back_edge_idx).set_class(weight);

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

  // TODO: what is b_it
  std::list<Bracket>::iterator b_it =
	this->back_edges.at(backedge_idx).bracket_it();
	// BracketList& bl = this->bracket_lists.at(vertex);

  // auto res = this->bracket_lists.at(vertex).erase(b_it);

  for (auto br = this->bracket_lists.at(vertex).begin();
	   br != this->bracket_lists.at(vertex).end();
	   br++) {
	if (br->back_edge_id() == b.id()) {
	  this->bracket_lists.at(vertex).erase(br);
	  //std::size_t s =  this->bracket_lists.at(vertex).size();
	  //this->bracket_lists.at(vertex).front().set_recent_size(s);
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
	Bracket(be.id(),
			core::constants::UNDEFINED_SIZE_T,
			core::constants::UNDEFINED_SIZE_T,
			be.is_capping_backedge());
  this->bracket_lists.at(vertex).push_front(br);
  std::list<Bracket>::iterator br_it = this->bracket_lists.at(vertex).begin();

  // std::size_t s = this->bracket_lists.at(vertex).size();
  // this->bracket_lists.at(vertex).set_size(++s);

  // update the backedge
  be.set_bracket_it(br_it);
  // TODO: unnecssary?
  //be.set_class(this->equiv_class_count_);
  //be.set_recent_size(s);
}

Bracket& Tree::top(std::size_t vertex) {
  return this->bracket_lists.at(vertex).front();
}

std::size_t Tree::new_class() { return this->equiv_class_count_++; }

Edge& Tree::get_incoming_edge(std::size_t vertex) {
  std::size_t e_idx = this->nodes.at(vertex).get_parent_idx();
  return this->tree_edges.at(e_idx);
}

/**
 * Get the node id of the node with the given sort value
 *
 * The sort vector is sorted in ascending order based on the index from 0..n
 * the value at index zero in the sort vector (`sort_`) is the node id of the
 * node with the smallest sort value and so forth.
 * We can then use this node id to get the node from the nodes vector (`nodes`)
 *
 * @param[in] idx the sort value of the node
 * @return the node id of the node with the given sort value
 */
std::size_t Tree::get_sorted(std::size_t idx) { return  this->sort_.at(idx);}

std::size_t Tree::get_sorted_g(std::size_t idx) { return  this->sort_g.at(idx);}
  
void Tree::cycles_vector(
  std::vector<std::tuple< size_t , size_t, size_t>>& vertices_four,
  std::vector<size_t> & classes
  ) {


  // a lambda that given a vertex id v will loop through the sort_ vector
  // and find the value of v in the sort_ vector and return the index
  // of that value
  auto find_sort_value = [this](std::size_t v) {
		for (std::size_t i = 0; i < this->sort_.size(); i++) {
			if (this->sort_.at(i) == v) {
				return i;
			}
		}
		return core::constants::UNDEFINED_SIZE_T;
  };

  // set of eq classes
  // if they have been seen or not
  std::set<size_t> seen;

  // number of vertices
  // the index is i & the value is the first/last index that i was seen
  std::vector<size_t> first_seen(this->size(), core::constants::UNDEFINED_SIZE_T);
  std::vector<size_t> last_seen(this->size(), core::constants::UNDEFINED_SIZE_T);

  // a stack of branching paths
  // contains node ids of nodes with more than 1 child
  std::stack<size_t> branching_stack; // TODO: rename branching paths?

  // a stack of suspended classes or vertices
  std::stack<size_t> suspended_stack;

  // which branches have been explored?
  // a map of vertex ids to a set of eq classes
  // the eq classes are the ones that have been explored
  std::map<size_t, std::set<size_t>> explored_branches;


  // a vector of class ids in order
  // std::vector<size_t> classes;
  // a vector of vertex ids in order
  std::vector<size_t> vertices;


  // a vertex of size_t pairs
  //std::vector<std::tuple< size_t , size_t, size_t>> vertices_four;


  /*
	initialize variables
   */

  // initialize vertices with id 0
  vertices.push_back(0);

  // current index to i is 1
  std::size_t i = 1;

  // the vertex at index i
  std::size_t v{core::constants::UNDEFINED_SIZE_T};

  // Is this a suspended run?
  bool suspended_run; // TODO: rename to suspended?

  // has i been set
  bool i_set{false};

  // initialize the current class
  std::size_t current_class{core::constants::UNDEFINED_SIZE_T};

  // outgoing backedge indexes
  std::set<size_t> obe_idxs;

  std::vector<Edge> children;

  // a lambda that will update classes and vertices with the current class and v values
  auto update_classes_and_vertices = [&]() {
	classes.push_back(current_class);
	vertices.push_back(v);


	std::size_t src = this->get_parent_edge(v).get_parent();
	std::size_t tgt = this->get_parent_edge(v).get_child();

	// push src and tgt to vertices_four
	vertices_four.push_back(std::make_tuple(current_class, src, tgt));


	// populate vertices_three

	// create a tuple of v and a set of classes for the current vertex and push it to vertices_three

	// loop through the classes in children and accumulate them into a set
	std::set<size_t> classes_set;
	for (Edge e : children) {
	  classes_set.insert(e.get_class());
	}



	// loop through the outgoing backedges of v and add the classes of these
	//  back-edges to the classes vector
	for (size_t obe_idx : obe_idxs) {

	  // ignore capping backedges
	  BackEdge obe = this->get_backedge(obe_idx);
	  if (!obe.is_capping_backedge()) {
		classes.push_back(obe.get_class());

		src = obe.get_src();
		tgt = obe.get_tgt();

		vertices_four.push_back(std::make_tuple(obe.get_class(), src, tgt));
	  }
	}
  };

  while (true) {

	suspended_run = false;
	v = this->get_sorted(i);
	current_class = this->get_parent_edge(v).get_class();

	// if the current class has been seen
	if (seen.count(current_class)) {
	  if (last_seen.at(current_class) == i) { suspended_run = true; }
	  last_seen.at(current_class) = i;
	}
	// it is the first time seeing this class
	else {
	  first_seen.at(current_class) = i;
	  seen.insert(current_class);
	}

	// get the vector of children
	children = this->get_child_edges(v);

	// is branching
	if (children.size() > 1) { branching_stack.push(v); }

	i_set = false;

	// obe_idxs

	obe_idxs = this->get_obe_idxs(v);

	// if node has outgoing back-edges and this is not the first time the class is being seen

	//std::cout << "suspended run: " << suspended_run
	//		  << " " << (first_seen[current_class] != i)
	//		  << " " << obe_idxs.empty() << std::endl;

	if (!suspended_run
		&& first_seen[current_class] != i
		&& !obe_idxs.empty()
	  ) {
	  for (size_t obe_idx : obe_idxs) {
		//std::cout << "obe_idx: " << obe_idx << std::endl;
		BackEdge obe = this->get_backedge(obe_idx);
		//obe.is_capping_backedge()
		// if outgoing backedge goes above the first seen or to zero then everything
		// (the subtree) below first seen is within the SESE bubble
		if (!obe.is_capping_backedge()
			&& obe.get_tgt() < first_seen[current_class]
			&& !branching_stack.empty()
		  ) {

		  std::size_t b = branching_stack.top();
		  // loop through the children of the nodes at b
		  for (Edge e : this->get_child_edges(b)) {
			// the child lies in v vertex id and not in i sort id
			std::size_t child = e.get_child();

			// if the child is not in explored branches then add it and then
			// push i into suspended_stack
			// TODO: does the lack of order break things here?
			if (!explored_branches[b].count(child)) {
			  explored_branches[b].insert(child);
			  suspended_stack.push(i);
			  i = find_sort_value(child);
			  i_set = true; // assumes find_sort_value above will always find a value
			  break;
			}
		  }
		}
	}
  }

	// if i is set continue
	if (i_set) { continue; }
	// if i is the last vertex
	else if (i == this->size() - 1) {
	  update_classes_and_vertices();
	  break;
	}
	// is branching
	else if (children.size() > 1) {
	  update_classes_and_vertices();
	  //	  std::cout << "here" << std::endl;
	  // loop through the edges in children and if you find one that is not in explored branches
	  // then push set i to that child and insert it into explored branches
	  for (Edge e : children) {
		std::size_t child = e.get_child();
		//std::cout << "child: " << child << std::endl;
		if (!explored_branches[v].count(child)) {
		  explored_branches[v].insert(child);
		  i = find_sort_value(child);
		  break;
		}
	  }
	}
	//
	else if (first_seen[current_class] == i // we are seeing the class for the first time
			 && !children.empty() // is not a leaf
	  ) {
	  update_classes_and_vertices();

	  std::size_t child = children.front().get_child();
	  i = find_sort_value(child);
	}
	// is a leaf or suspended run is true
	else if (children.empty() || suspended_run) {
	  update_classes_and_vertices();


	  std::size_t b = branching_stack.top();

	  // loop through the explored branches of v and if all have been explored then pop
	  // the branching stack then set i to the top of the suspended stack and pop the suspended stack
	  // if the suspended stack is empty then we are done
	  if (explored_branches[b].size() == this->get_child_edges(b).size()) {
		branching_stack.pop();
		// TODO: what to do in this case?
		if (suspended_stack.empty()) {
		  std::cout << "suspended stack is empty, exiting" << std::endl;
		  break;
		}
		i = suspended_stack.top();
		suspended_stack.pop();
	  }
	  else {
		// take the value on top of the branching stack and find the first child that has not been explored
		std::size_t b = branching_stack.top();
		for (Edge e : this->get_child_edges(b)) {
		  std::size_t child = e.get_child();
		  if (!explored_branches[b].count(child)) {
			explored_branches[b].insert(child);
			i = find_sort_value(child);
			//if (!suspended_run) {  }
			break;
		  }
		}

		if (explored_branches[v].size() == children.size()) {
		  branching_stack.pop();
		  //if (suspended_run) { i = suspended_stack.top(); }
		  //if (!suspended_stack.empty()) { suspended_stack.pop(); }

		}
	  }
	}
	// i is not branching
	// same as 4
	else if (children.size() == 1) {
	  update_classes_and_vertices();

	  std::size_t child = children.front().get_child();
	  i = find_sort_value(child);
	}
	// suspended stack is empty
	else if (suspended_stack.empty()) {
	  // loop through the edges in children and if you find one that is not in explored branches
	  // then push set i to that child and insert it into explored branches
	  for (Edge e : children) {
		std::size_t child = e.get_child();
		if (!explored_branches[v].count(child)) {
		  explored_branches[v].insert(child);
		  i = find_sort_value(child);
		  break;
		}
	  }
	}
	// if i is not set
	else if (!i_set) {
	   // suspended is assumed to not be empty
	   // happens if we have explored all paths
	  branching_stack.pop();
	  i = suspended_stack.top();
	  //if (suspended_run) {  }
	  suspended_stack.pop();
	}
	else {
	  std::cout << "something went wrong" << std::endl;
	  break;
	}

  } // end while loop


  if (false) {
	// print class vector
	std::cout << "classes: ";
	for (size_t c : classes) {
	  std::cout << c << " ";
	}
	std::cout << std::endl;

	// print vertices_four vector of pairs
	std::cout << std::endl << "vertices_four: ";
	for (auto v : vertices_four) {
	  std::cout << "(" << std::get<0>(v) << ": " << std::get<1>(v)
				<< ", " << std::get<2>(v) << "), ";
	}

	std::cout << std::endl;
  }
}

std::vector<Edge> Tree::compute_edge_stack() {
  std::vector<Edge> edge_stack;
  std::size_t current_vertex{};
  std::set<size_t> seen;
  std::set<size_t> exp;

  std::stack<std::size_t> bts;

  int counter {};

  while (counter < 15) {

	std::cout << "current vertex: " << current_vertex << std::endl;

	++counter;
	seen.insert(current_vertex);

	std::size_t prev_vertex = current_vertex;


	std::vector<Edge> child_edges = this->get_child_edges(current_vertex);


	std::set<size_t> o = this->get_obe(current_vertex);

	if (!child_edges.empty() && !o.empty()) {
	  bts.push(current_vertex);
	}

	for (auto it : child_edges) {
	  auto c = it.get_child();

	  if (!seen.count(c)) {
		current_vertex = c;
		edge_stack.push_back(it);
		exp.insert(c);
		break;
	  }

	}

	for (auto it : child_edges) {
	  auto c = it.get_child();



	  if (!exp.count(c)) {
		bts.push(prev_vertex);
		break;
	  }
	}

	// print child edges for current vertex
	//for (auto e : child_edges) {
	// std::cout << e.get_parent() << " " << e.get_class() << " " << e.get_child() << std::endl;
	  //}



	if (child_edges.empty()) {
	  current_vertex = bts.top();
	  bts.pop();
	}
  }

  std::cout << "edge stack" << std::endl;

  for (auto e : edge_stack ) {
	std::cout << e.get_parent() << " " << e.get_class() << " " << e.get_child() << std::endl;
  }

  return edge_stack;
}

std::vector<std::pair<std::size_t, std::size_t>>
  Tree::compute_edge_stack2() {

  std::stack<std::pair<std::size_t, std::size_t>> s;
  // first is vertex, second is eq class of the vertex
  std::vector<std::pair<std::size_t, std::size_t>> v;
  
  for (std::size_t j{}; j < this->size(); j++) {
	std::size_t i = this->get_sorted(j);
	std::size_t k = this->get_sorted_g(i);

	// the 0th vertex has no parent and will throw an exception
	// if we try to get its parent
	// this should be fixed somehow but for now we just skip it
	if (k > 0) {
	  	Edge const& parent_edge = this->get_parent_edge(k);
		//std::cout << "parent: " << parent_edge.get_class()
		//		  << " " << parent_edge.get_color()
		//		  << std::endl;
		if (parent_edge.get_color() == core::color::black) {
		  v.push_back(std::make_pair(k, parent_edge.get_class()));  
		  s.push(std::make_pair(k, parent_edge.get_class()));  
		}
		
	}

	std::set<size_t> obes = this->get_obe_idxs(k);
	for (auto o : obes) {
	  BackEdge& be  = this->get_backedge(o);
	  if (be.is_capping_backedge()) { continue; }	  		
	  //std::cout << "backedge: " << be.get_class()
	  //			<< " " << be.get_color()
	  //			<< std::endl;

	  if (be.get_color() == core::color::black) {
		v.push_back(std::make_pair(k, be.get_class()));  
		  s.push(std::make_pair(k, be.get_class()));  
		}
	}		
  }

  // print the conents of the stack
  for (auto it : v) {
	std::cout << it.first << " " << it.second << std::endl;
	}

  return v;
}
  
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

	std::size_t i = this->get_sorted(j);

	// Vertex const& vv = this->get_vertex(i);
	for (auto c : this->get_child_edges(i)) {
	  std::string cl = c.get_class_idx() > 10000
						   ? "\u2205"
						   : std::to_string(c.get_class_idx());

	  std::string color = c.get_color() == core::color::gray ? "gray" : "black";

	  std::cout << std::format("\t{}  -- {}  [label=\"{} {}\" color=\"{}\"];\n",
							   i, c.get_child(), c.id(), cl, color);
	}

	//for (auto c : this->get_children_w_id(i)) {
	  //std::cout << std::format("\t{} -- {}  [label=\"{}\"];\n", i, c.second, c.first);
	//}

	for (auto o: this->get_obe_w_id(i)) {
	  std::string cl = o.first > 10000 ? "\u2205" : std::to_string(o.first);
	  // a capping backedge is red and can never have been gray
	  BackEdge be = this->get_backedge_given_id(o.first);

	  //std::cout << "be: " << be.get_src() << " " << be.get_tgt() << " " << o.first << " "  << o.second << "\n";

	  bool is_capping = be.is_capping_backedge();

	  std::string class_ = be.get_class() ==  core::constants::UNDEFINED_SIZE_T
		 ?  "" : std::to_string(be.get_class());

	  std::string color{};

	  if (is_capping) {
		color = "red";
	  } else if (this->get_backedge_given_id(o.first).get_color() == core::color::gray) {
		color = "gray";
	  } else if (this->get_backedge_given_id(o.first).get_color() == core::color::black) {
		color = "black";
	  } else {
		color = "blue";
	  }

	  std::cout << std::format(
		  "\t{} -- {} [label=\"{} {}\" style=\"dotted\" penwidth=\"3\" color=\"{}\"];\n",
		  i, be.get_tgt(), cl, class_, color
		  );
	}
  }

  std::cout << "}" << std::endl;
}

} // namespace spanning_tree
