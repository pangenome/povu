#ifndef SPANNING_TREE_HPP
#define SPANNING_TREE_HPP

#include <cstddef>
#include <string>
#include <vector>
#include <set>
#include <list>
#include <map>
#include <unordered_map>

#include "../common/types.hpp"
#include "./bracket_list.hpp"

namespace povu::spanning_tree {
using namespace povu::graph_types;
using namespace povu::bracket_list;
namespace pgt = povu::graph_types;


// prototype the classes
class Edge;
class Vertex;
class BackEdge;

enum class EdgeType {
  tree_edge,
  back_edge,
  capping_back_edge,
  simplifying_back_edge,
};

/*
 * edge
 * ----
 *
 * a tree edge
 */
class Edge {
  std::size_t id_; // (recent class)

  // rename to parent and child
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
  // Bracket* b; // if it is a backedge
  bool null_;

  color color_;

public:
  // --------------
  // constructor(s)
  // --------------
  Edge();
  Edge(std::size_t id, std::size_t src, std::size_t tgt, color c=color::black);

  // ---------
  // getter(s)
  // ---------
  std::size_t id() const;
  std::size_t get_child() const; // TODO: remove, use get_child_v_idx
    /**
    * @brief get the index of the parent vertex
    *
    * @return the index of the parent vertex
   */
  std::size_t get_child_v_idx() const;
  color get_color() const;
  std::size_t get_parent() const; // TODO: remove, superceded by get_parent_v_idx
  /**
    * @brief get the index of the parent vertex
    *
    * @return the index of the parent vertex
   */
  std::size_t get_parent_v_idx() const;

  std::size_t get_v1() const;
  std::size_t get_v2() const;

  // FIXME: we need both because of a non const call that depends on the non const
  std::size_t get_class() const;
  std::size_t get_class_idx();

  // ---------
  // setter(s)
  // ---------
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


  // TODO: why this duplication?
  std::size_t class_; // equivalnce class id
  std::size_t recent_class_; // FIXME: remove this or use it
  std::size_t recent_size_; //

  // TODO: use an enum here
  //bool capping_back_edge_; // is a capping back edge // TODO: remove, use type_
  EdgeType type_;
  bool null_; // FIXME: remove this or use it

  color color_;

public:
  BackEdge(std::size_t id, std::size_t src, std::size_t tgt, EdgeType t, color c=color::black);

  std::size_t id() const;
  std::size_t get_src() const;
  std::size_t get_tgt() const;

  std::size_t get_class() const;
  bool is_class_defined() const;

  bool is_capping_backedge() const;
  color get_color() const;

  void set_class(std::size_t c);

  EdgeType type() const;
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

  // TODO: rename to something better?
  std::size_t g_v_id_ {}; // id/name of the vertex in the input GFA

  pgt::VertexType type_;

  /*
   dfs_num of the highest node originating from an outgoing backedge from this
   vertex or from a child of this vertex
   */
  std::size_t hi_;

  bool null_;

public:
  // --------------
  // constructor(s)
  // --------------
  //Vertex(); // creates a root with parent id set to and id of zero
  //Vertex(std::size_t id, std::size_t parent_id);
  Vertex(std::size_t dfs_num, std::size_t g_v_id, VertexType type_);
  // ---------
  // getter(s)
  // ---------
  bool is_root() const;
  bool is_leaf() const;
  std::size_t dfs_num() const;
  std::size_t parent() const; // TODO: remove
  std::size_t hi() const; // TODO: remove
  std::size_t g_v_id() const;
  VertexType type() const;

  std::set<size_t> const& get_obe() const;
  std::set<size_t> const& get_ibe() const;
  bool is_null() const;

  // get the index of the edge that points to the parent in the tree

  size_t const& get_parent_idx() const; // TODO: remove, superceded by get_parent_edge_idx
  size_t get_parent_e_idx() const;

  std::set<size_t> const& get_children() const;

  // ---------
  // setter(s)
  // ---------
  void add_obe(std::size_t obe_id);
  void add_ibe(std::size_t ibe_id);
  void add_child(std::size_t e_id);
  void unset_null(); // sets null to false;

  // the index of the parent node in the tree vertex
  void set_parent(std::size_t n_id);
  void set_g_v_id(std::size_t g_v_id);
  void set_type(VertexType t);
  void set_hi(std::size_t val);
  // the dfs num of the node
  void set_dfs_num(std::size_t idx);
};

class Tree {
  // no of nodes in the tree
  std::vector<Vertex> nodes;
  std::vector<Edge> tree_edges;
  std::vector<BackEdge> back_edges;

  // TODO: move the bracket list don't keep creating a new one
  // a BracketList for each node
  // the list of backedges bracketing a node
  // a bracketList is a backedge with some metadata around it
  //std::vector<BracketList *> bracket_lists;
  std::vector<WBracketList*> bracket_lists;

  // the index of an edge or backedge in the input graph
  // key is the graph_idx and the value is the tree_idx or the backedge idx
  std::map<std::size_t, std::pair<EdgeType, std::size_t>> g_edge_idx_map;

  // tree idx to graph edge idx
  std::map<std::size_t, std::size_t> tree_graph_idx_map_;

  // a map from the back edge id to the index in the back_edges vector
  std::map<std::size_t, std::size_t> be_id_to_idx_map_;

  // a map from the back edge id to the index in the back_edges vector
  std::map<std::size_t, std::size_t> te_id_to_idx_map_;

  /*
    a map from the edge id (in the spanning tree) to the index in the tree_edges
    vector or the back_edges vector
   */
  std::map<std::size_t, std::pair<EdgeType, std::size_t>> edge_id_map_;

  // Holds the topo mapping of the tree
  // the index is the position in the toposort
  // the value at a poistion is the index in nodes
  // topo sort vector
  std::vector<std::size_t> sort_;

  // sort based on the input graph
  // the index in the vertex is the index in the input graph
  // and the value is the index in the tree
  std::vector<std::size_t> sort_g;

  // TODO: replace sort and sort_g with a two way map or remove both

  static const size_t root_node_index {}; // 0


  std::size_t equiv_class_count_;

public:

  // --------------
  // constructor(s)
  // --------------
  Tree(std::size_t size);

  // ---------
  // destructor(s)
  // ---------
  ~Tree();

  // ---------
  // getter(s)
  // ---------
  Vertex& get_root();
  std::size_t get_root_idx() const;

  // number of vertices in the tree
  std::size_t size() const;
  std::size_t tree_edge_count() const;
  std::size_t back_edge_count() const;


  Vertex const &get_vertex(std::size_t vertex) const;
  Vertex& get_vertex_mut(std::size_t vertex);
  // get the parent vertex of the vertex at v_idx
  // if v_idx is the root, caues undefined behaviour
  Vertex const& get_p_vtx(std::size_t v_idx) const;


  // given  tree edge index return a the child
  // returns a node index
  // a set?


  // takes a vertex index and returns pair:
  //   first => backedge id
  //   second =>  node idx of the target vertex
  std::set<std::pair<std::size_t, std::size_t>> get_obe_w_id(std::size_t vertex);
  std::set<std::pair<std::size_t, std::size_t>> get_ibe_w_id(std::size_t vertex);

  // take a vertex index and return an edge id and the child node index
  std::set<std::pair<std::size_t, std::size_t>> get_children_w_id(std::size_t vertex);


  // return edges that point to children of the vertex
  std::vector<Edge> get_child_edges(std::size_t vertex);


  // get index of the  be in back_edges vector
  std::set<std::size_t> get_obe_idxs(std::size_t vertex);
  std::set<std::size_t> get_ibe_idxs(std::size_t vertex);

  size_t list_size(std::size_t vertex);
  size_t get_hi(std::size_t vertex);


  /**
    * @brief get indexes of the vertices the obes from this vertex points to (tgt/targets)
   */
    // TODO rename to get_obe_tgts
  std::set<size_t> get_obe(std::size_t vertex); // get backedge target indexes

  /**
    * @brief get sources of the vertices the ibes from this vertex points from (srcs)
   */
      // TODO rename to get_obe_srcs
  std::set<size_t> get_ibe(std::size_t vertex);

  /**
   * @brief a reference to the tree edge given the index in the tree_edges vector
   *
   * @param edge_idx the index of the edge in the tree_edges vector
   * @return a reference to the edge
   */
  const Edge& get_tree_edge(std::size_t edge_idx) const;

  /**
    * @brief a reference to the tree edge given the edge id
    *
    * @param edge_id the id of the edge
    * @return a reference to the edge
   */
  //const Edge& get_tree_edge_by_id(std::size_t edge_id) const;

  std::size_t get_graph_edge_id(std::size_t tree_edge_id) const;

  const std::pair<EdgeType, std::size_t>& get_edge_idx(std::size_t edge_id) const;

  // return reference to a back edge given the
  // index of the back edge in the back_edges vector
  BackEdge& get_backedge(std::size_t backedge_idx);
  BackEdge& get_backedge_ref_given_id(std::size_t backedge_id);
  // given the back edge's unique back edge id return a reference to the backedge
  BackEdge get_backedge_given_id(std::size_t backedge_id);
  std::set<std::size_t> get_children(std::size_t vertex);

  bool is_root(std::size_t vertex) const;
  bool is_leaf(std::size_t vertex) const;

  // given a vertex id, return a reference to the edge that points to the parent
  Edge const& get_parent_edge(std::size_t vertex) const;


  // TODO: rename to get_parent_edge_mut
  Edge& get_incoming_edge(std::size_t vertex);

  /**
    * @brief get the index of the edge to the parent vertex
    * if v_idx is the root index it causes to undefined behavior
   */
  std::size_t get_parent(std::size_t v_idx); // TODO: remove, superceded by get_parent_idx
  std::size_t get_parent_v_idx(std::size_t v_idx) const;

  bool has_ibe(std::size_t vertex, std::size_t qry_idx);
  bool has_obe(std::size_t vertex, std::size_t qry_idx);
  bool has_child(std::size_t vertex, std::size_t qry_idx);

  // the vertex id of the node at sort value idx
  std::size_t get_sorted(std::size_t idx);

  std::size_t get_sorted_g(std::size_t idx);

  const std::map<std::size_t, std::pair<EdgeType, std::size_t>>& get_g_edge_idx_map() const;

  // -------
  // setters
  // -------

  //void topo_sort();

  // takes an index in toposort
  // and a vertex and sets the vertex as value in the toposort
  // vector
  void set_sort(std::size_t idx, std::size_t vertex);

  // takes an index in the input graph
  // and a vertex and sets the vertex as value in the toposort
  // vector
  void set_sort_g(std::size_t idx, std::size_t vertex);


  void add_vertex(Vertex&& v);

  // set the dfs number of a vertex
  void set_dfs_num(std::size_t vertex, std::size_t dfs_num);
  void set_vertex_type(std::size_t vertex, VertexType type);

  std::size_t add_be(std::size_t frm, std::size_t to, EdgeType t, color clr=color::black);

  std::size_t add_be(std::size_t frm,
                     std::size_t to,
                     std::size_t g_edge_id,
                     EdgeType t,
                     color clr=color::black);

  void add_tree_edge(std::size_t frm,
                     std::size_t to,
                     std::size_t g_edge_idx,
                     color clr=color::black);

  void set_hi(std::size_t vertex, std::size_t val);

  // vst
  void concat_bracket_lists(std::size_t parent_vertex, std::size_t child_vertex);
  void del_bracket(std::size_t vertex, std::size_t backedge_idx);
  void push(std::size_t vertex, std::size_t backege_idx);

  // get the bracket on top of the bracket list for v
  Bracket& top(std::size_t vertex);

  // return the current equivalence class count then increment it
  std::size_t new_class();

  BracketList& get_bracket_list(std::size_t vertex);

  // ------------
  // I/O
  // ------------

  void print_dot();
};

} // namespace spanning_tree


#endif
