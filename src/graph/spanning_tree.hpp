#ifndef SPANNING_TREE_HPP
#define SPANNING_TREE_HPP

#include <cstddef>
#include <vector>
#include <set>
#include <memory>

namespace spanning_tree {
// prototype the classes
class Edge; class Vertex; class Bracket; class BracketList;

// enum colour { grey, black };

/*
 * spanning tree
 * -------------
 */
class Edge {
  std::size_t id_; // (recent class)

  std::size_t src; // target vertex
  std::size_t tgt; // source vertex

  // does a tree edge need color?
  // colour edge_colour;

  std::size_t class_; // equivalnce class id

  // TODO: is size used?
  std::size_t size;      // (recent size) size of the bracket list
  // std::size_t recent_class; // (recent class) id of the topmost backedge

  std::size_t backedge_id; // TODO: used??
  // std::size_t id; // (recent class)
  //Bracket* b; // if it is a backedge
  bool null_;
public:

  Edge();
  Edge(std::size_t id, std::size_t src, std::size_t tgt);


  std::size_t id() const;
  std::size_t get_parent() const;
  std::size_t get_child() const;

  // TODO: remove _idx
  std::size_t get_class_idx();
  void set_class_idx(std::size_t c);
};

class BackEdge {
  std::size_t id_; //
  std::size_t src; // target vertex
  std::size_t tgt; // source vertex

  Bracket* b; // pointer to bracket in the BracketList
  std::size_t class_; // equivalnce class id
  std::size_t recent_class_; //
  std::size_t recent_size_; //

  bool capping_back_edge_; // is a capping back edge
  bool null_;

public:
  BackEdge(bool capping_be=false); // TODO remove? is this used?
  BackEdge(std::size_t id, std::size_t src, std::size_t tgt, bool capping_be=false);

  std::size_t id() const;
  std::size_t get_src() const;
  std::size_t get_tgt() const;

  std::size_t get_class() const;
  bool is_class_defined() const;
  std::size_t get_recent_class() const;
  std::size_t get_recent_size() const;

  bool is_capping_backedge() const;

  void set_class(std::size_t c);
  void set_recent_class(std::size_t c);
  void set_recent_size(std::size_t s);

  Bracket* bracket_ptr();
  void set_bracket_ptr(Bracket* p);
};

/*
 * Bracket
 * -------
 */
class Bracket {
  std::size_t id; // backedge id? TODO: remove?
  std::size_t recent_size_;
  std::size_t recent_class_; // TODO: rename to class?
  // std::size_t back_edge_id_;
  Bracket *next_;
  Bracket *prev_;

public:
  Bracket();
  Bracket(BackEdge& b, std::size_t recent_size, std::size_t recent_class);

  Bracket* next();
  Bracket* prev();
  std::size_t recent_size() const;
  std::size_t recent_class() const;

  void set_next(Bracket *p);
  void set_prev(Bracket *p);
  void set_recent_size(std::size_t s);
  void set_recent_class(std::size_t c);
};

/*
 * Bracket List
 * ------------
 *
 * Maintains a list of brackets
 */
class BracketList {
  Bracket *first_;
  Bracket *last_;

  std::size_t node_id_;
  std::size_t size_;

public:
  BracketList();

  std::size_t size() const;

  Bracket& top();
  void push(Bracket* e);

  void set_size(std::size_t s);

  Bracket* first();
  Bracket* last();

  // delete
  void remove();

  // concat
  //void append_ptr(BracketList* b_ptr);
  void prepend(BracketList& b);
  void append(BracketList& b);
};

// TODO: remove unused methods
class Vertex {
  std::size_t idx; // dfsnum an id in toposort. Is this necessary?

  std::size_t parent_id; // id to idx // index to the tree edge vector ?
  // indexes of the children edges in the tree_edges vector
  std::set<std::size_t> children; // children // index to the tree edge vector

  std::set<size_t> obe; // out back edges
  std::set<size_t> ibe;  // in back edges

  /*
   id of the highest node originating from an outgoing backedge from this
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
  std::size_t id() const;
  std::size_t parent() const; // TODO: remove
  std::size_t hi() const; // TODO: remove

  std::set<size_t> const& get_obe() const;
  std::set<size_t> const& get_ibe() const;
  // get the index of the edge that points to the parent in the tree
  size_t const& get_parent_idx() const; // TODO: get_parent_edge_idx
  std::set<size_t> const& get_children() const;

  void add_obe(std::size_t id);
  void add_ibe(std::size_t id);
  void add_child(std::size_t id);
  void unset_null(); // sets null to false;
  void set_parent(std::size_t id);
  void set_hi(std::size_t val);
};

class Tree {
  std::vector<Vertex> nodes;
  std::vector<Edge> tree_edges;
  std::vector<BackEdge> back_edges;
  std::vector<BracketList> bracket_lists;

  static const size_t root_node_index{0};

  // std::size_t list_size_; // TODO: remove
  std::size_t equiv_class_count_;

public:
  // constructor(s)
  Tree();
  Tree(std::size_t size);

  // getters
  Vertex& get_root();
  std::size_t size() const;

  Vertex const &get_vertex(std::size_t vertex) const;

  Edge& get_incoming_edge(std::size_t vertex);

  // given  tree edge index return a the child
  // returns a node index
  // a set?
  std::set<std::pair<std::size_t, std::size_t>> get_obe_w_id(std::size_t vertex);
  std::set<std::pair<std::size_t, std::size_t>> get_ibe_w_id(std::size_t vertex);
  std::set<std::pair<std::size_t, std::size_t>> get_children_w_id(std::size_t vertex);

  std::vector<Edge> get_child_edges(std::size_t vertex);
  std::set<std::size_t> get_obe_idxs(std::size_t vertex); // get index of obe in back_edges vector
  std::set<std::size_t> get_ibe_idxs(std::size_t vertex); // get index of obe in back_edges vector

  size_t list_size(std::size_t vertex);
  size_t get_hi(std::size_t vertex);
  std::set<size_t> get_obe(std::size_t vertex); // get backedge target indexes
  std::set<size_t> get_ibe(std::size_t vertex);
  std::set<std::size_t> get_children(std::size_t vertex);

  bool is_root(std::size_t vertex) const;
  std::size_t get_parent(std::size_t vertex);
  BackEdge& get_backedge(std::size_t backedge_idx);
  bool has_ibe(std::size_t vertex, std::size_t qry_idx);
  bool has_obe(std::size_t vertex, std::size_t qry_idx);
  bool has_child(std::size_t vertex, std::size_t qry_idx);

  // setters
  void add_vertex(Vertex&& v);
  std::size_t add_be(std::size_t frm, std::size_t to, bool capping_be=false);
  void add_tree_edge(std::size_t frm, std::size_t to);

  void set_hi(std::size_t vertex, std::size_t val);

  // vst
  void concat_brackets(std::size_t parent_vertex, std::size_t child_vertex);
  void del_bracket(std::size_t vertex, std::size_t backedge_idx);
  void push(std::size_t vertex, std::size_t backege_idx);

  // get the bracket on top of the bracket list for v
  Bracket& top(std::size_t vertex);

  // return the current equivalance class count then increment it
  std::size_t new_class();

  // I/O
  void print_dot();
};

} // namespace spanning_tree
#endif
