#ifndef SPANNING_TREE_HPP
#define SPANNING_TREE_HPP

#include <cstddef>
#include <vector>
#include <set>
#include <memory>
#include <list>


#include "../core/core.hpp"

namespace spanning_tree {
// prototype the classes
class Edge;
class Vertex;
class Bracket; // holds metadata about a back edge
class BackEdge; //
typedef std::list<Bracket> BracketList;


/*
 * edge
 * ----
 *
 * a tree edge
 */
class Edge {
  std::size_t id_; // (recent class)

  std::size_t src; // target vertex
  std::size_t tgt; // source vertex

  // does a tree edge need color?
  // color edge_color;

  std::size_t class_; // equivalnce class id

  // TODO: is size used?
// not used in tree edge
  std::size_t size;      // (recent size) size of the bracket list
  // std::size_t recent_class; // (recent class) id of the topmost backedge

  std::size_t backedge_id; // TODO: used??
// not used in tree edge
  // std::size_t id; // (recent class)
  //Bracket* b; // if it is a backedge
  bool null_;

  core::color color_;

public:

  Edge();
  Edge(std::size_t id, std::size_t src, std::size_t tgt, core::color c=core::color::black);

  std::size_t id() const;
  std::size_t get_parent() const;
  std::size_t get_child() const;

  core::color get_color() const;

  // FIXME: we need both because of a non const call that depends on the non const
  std::size_t get_class() const;
  std::size_t get_class_idx();

  void set_class_idx(std::size_t c); // deprecated replaced by set_class
  void set_class(std::size_t c);
};

/*
 * back edge
 * ----------
 *
 * an edge from a vertex to an ancestor (not parent) in the spanning tree
 */
class BackEdge {
  std::size_t id_; // a unique indeifier of the backedge
  std::size_t src; // target vertex
  std::size_t tgt; // source vertex

  // iterator to bracket in the bracketlist
  std::list<Bracket>::iterator bi;

  //Bracket* b; // pointer to bracket in the BracketList
  std::size_t class_; // equivalnce class id
  std::size_t recent_class_; //
  std::size_t recent_size_; //

  bool capping_back_edge_; // is a capping back edge
  bool null_;

  core::color color_;

public:
  BackEdge(bool capping_be=false); // TODO remove? is this used?
  BackEdge(std::size_t id, std::size_t src, std::size_t tgt, bool capping_be=false, core::color c=core::color::black);

  std::size_t id() const;
  std::size_t get_src() const;
  std::size_t get_tgt() const;

  std::size_t get_class() const;
  bool is_class_defined() const;
  std::size_t get_recent_class() const;
  std::size_t get_recent_size() const;

  bool is_capping_backedge() const;
  core::color get_color() const;

  void set_class(std::size_t c);
  void set_recent_class(std::size_t c);
  void set_recent_size(std::size_t s);

  // TODO: check if it is the end?
  std::list<Bracket>::iterator bracket_it();
  void set_bracket_it(std::list<Bracket>::iterator it);
  // Bracket* bracket_ptr();
  // void set_bracket_ptr(Bracket* p);
};

/*
 * Bracket
 * -------
 *
 *
 */
class Bracket {
  std::size_t back_edge_id_; // rename to backedge id? TODO: remove?
  std::size_t recent_size_;
  std::size_t recent_class_; // TODO: rename to class?
  // std::size_t back_edge_id_;
  //Bracket *next_;
  //Bracket *prev_;

  //std::size_t class_; // equivalnce class id

  // if the brackedge this bracket represents is a capping backedge
  bool is_capping_;

public:
  Bracket();
  Bracket(std::size_t backedge_id, std::size_t recent_size, std::size_t recent_class, bool is_capping=false);

  std::size_t back_edge_id();

  // Bracket* next();
  // Bracket* prev();

  // true if the backedge described by this bracket is a capping backedge
  // else false
  bool is_capping() const;
  std::size_t recent_size() const;
  std::size_t recent_class() const;

  void set_recent_size(std::size_t s);
  void set_recent_class(std::size_t c);
  //void set_back_edge_class(std::size_t c);


  // void set_next(Bracket *p);
  // void set_prev(Bracket *p);

};

/*
 * Vertex
 * ------
 *
 * TODO: remove unused methods
 */
class Vertex {
  // dfsnum of the node in toposort
  std::size_t dfs_num_;

  std::size_t parent_id; // id to idx // index to the tree edge vector ?

  // indexes of the children edges in the tree_edges vector
  std::set<std::size_t> children; // children // index to the tree edge vector

  std::set<size_t> obe; // out back edges
  std::set<size_t> ibe;  // in back edges

  /*
   dfs_num of the highest node originating from an outgoing backedge from this
   vertex or from a child of this vertex
   */
  std::size_t hi_;

  bool null_;

public:
  // constructor(s)
  Vertex(); // creates a root with parent id set to and id of zero
  Vertex(std::size_t id, std::size_t parent_id);

  bool is_root() const;
  bool is_leaf() const;
  std::size_t dfs_num() const;
  std::size_t parent() const; // TODO: remove
  std::size_t hi() const; // TODO: remove

  std::set<size_t> const& get_obe() const;
  std::set<size_t> const& get_ibe() const;

  // get the index of the edge that points to the parent in the tree
  size_t const& get_parent_idx() const; // TODO: rename to get_parent_edge_idx

  std::set<size_t> const& get_children() const;

  void add_obe(std::size_t obe_id);
  void add_ibe(std::size_t ibe_id);
  void add_child(std::size_t e_id);
  void unset_null(); // sets null to false;

  // the index of the parent node in the tree vertex
  void set_parent(std::size_t n_id);

  void set_hi(std::size_t val);
  // the dfs num of the node
  void set_dfs_num(std::size_t idx);

};

class Tree {
  // no of nodes in the tree
  std::vector<Vertex> nodes;
  std::vector<Edge> tree_edges;
  std::vector<BackEdge> back_edges;

  // a BracketList for each node
  // the list of backedges bracketing a node
  // a bracketList is a backedge with some metadata around it
  std::vector<BracketList> bracket_lists;

  // Holds the topo mapping of the tree
  // the index is the position in the toposort
  // the value at a poistion is the index in nodes
  // topo sort vector
  std::vector<std::size_t> sort_;

  // sort based on the input graph
  // the index in the vertex is the index in the input graph
  // and the value is the index in the tree
  std::vector<std::size_t> sort_g;

  static const size_t root_node_index{0};

  // std::size_t list_size_; // TODO: remove
  std::size_t equiv_class_count_;

public:

  // --------------
  // constructor(s)
  // --------------
  Tree();
  Tree(std::size_t size);

  // ---------
  // getter(s)
  // ---------
  Vertex& get_root();
  //
  std::size_t size() const;
  Edge& get_incoming_edge(std::size_t vertex);

  Vertex const &get_vertex(std::size_t vertex) const;
  Vertex& get_vertex_mut(std::size_t vertex);


  // given  tree edge index return a the child
  // returns a node index
  // a set?


  // takes a vertex index and returns pair:
  //   first => backedge id
  //   second =>  node idx of the target vertex
  std::set<std::pair<std::size_t, std::size_t>> get_obe_w_id(std::size_t vertex);
  std::set<std::pair<std::size_t, std::size_t>> get_ibe_w_id(std::size_t vertex);

  // take a vertex index and return a
  std::set<std::pair<std::size_t, std::size_t>> get_children_w_id(std::size_t vertex);

  // return edges that point to children of the vertex
  std::vector<Edge> get_child_edges(std::size_t vertex);

  // given a vertex id, return a reference to the edge that points to the parent
  Edge const& get_parent_edge(std::size_t vertex) const;

  // get index of the  be in back_edges vector
  std::set<std::size_t> get_obe_idxs(std::size_t vertex);
  std::set<std::size_t> get_ibe_idxs(std::size_t vertex);

  size_t list_size(std::size_t vertex);
  size_t get_hi(std::size_t vertex);
  std::set<size_t> get_obe(std::size_t vertex); // get backedge target indexes
  std::set<size_t> get_ibe(std::size_t vertex);

  // return reference to a back edge given the
  // index of the back edge in the back_edges vector
  BackEdge& get_backedge(std::size_t backedge_idx);
  BackEdge& get_backedge_ref_given_id(std::size_t backedge_id);
  // given the back edge's unique back edge id return a reference to the backedge
  BackEdge get_backedge_given_id(std::size_t backedge_id);
  std::set<std::size_t> get_children(std::size_t vertex);

  bool is_root(std::size_t vertex) const;
  std::size_t get_parent(std::size_t vertex);

  bool has_ibe(std::size_t vertex, std::size_t qry_idx);
  bool has_obe(std::size_t vertex, std::size_t qry_idx);
  bool has_child(std::size_t vertex, std::size_t qry_idx);

  // the vertex id of the node at sort value idx
  std::size_t get_sorted(std::size_t idx);


  std::size_t get_sorted_g(std::size_t idx);

  // -------
  // setters
  // -------

  // takes an index in toposort
  // and a vertex and sets the vertex as value in the toposort
  // vector
  void set_sort(std::size_t idx, std::size_t vertex);

  // takes an index in the input graph
  // and a vertex and sets the vertex as value in the toposort
  // vector
  void set_sort_g(std::size_t idx, std::size_t vertex);

  
  // set the dfs number of a vertex
  void set_dfs_num(std::size_t vertex, std::size_t dfs_num);


  std::size_t add_be(std::size_t frm,
					 std::size_t to,
					 bool capping_be=false,
					 core::color color=core::color::black);

  std::size_t add_be(std::size_t frm,
					 std::size_t to,
					 std::size_t weight,
					 bool capping_be=false,
					 core::color color=core::color::black);

  void add_tree_edge(std::size_t frm,
					 std::size_t to,
					 core::color color=core::color::black);

  // returns the tree edge index
  std::size_t add_tree_edge(std::size_t frm,
							std::size_t to,
							std::size_t weight,
							core::color color=core::color::black);

  void set_hi(std::size_t vertex, std::size_t val);

  // vst
  void concat_bracket_lists(std::size_t parent_vertex, std::size_t child_vertex);
  void del_bracket(std::size_t vertex, std::size_t backedge_idx);
  void push(std::size_t vertex, std::size_t backege_idx);

  // get the bracket on top of the bracket list for v
  Bracket& top(std::size_t vertex);

  // return the current equivalance class count then increment it
  std::size_t new_class();

  BracketList& get_bracket_list(std::size_t vertex);

  // compute equiv classes stack and the vertices vector
  //void cycles_vector();

  void cycles_vector(std::vector<std::tuple< size_t , size_t, size_t>>& v, std::vector<std::size_t>& classes);

  std::vector<Edge> compute_edge_stack();

  // added after biedging
  std::vector<std::pair<std::size_t, std::size_t>> compute_edge_stack2();

  // I/O
  void print_dot();
};

} // namespace spanning_tree
#endif
