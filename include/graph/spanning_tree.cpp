#include "./spanning_tree.hpp"
#include <vector>

namespace povu::spanning_tree {

/*
 * Edge
 * ----
 */
/*  constructor(s) */
Edge::Edge(pt::id_t id, pt::idx_t parent_idx, pt::idx_t child_idx, pgt::color_e c)
  : id_(id), parent_idx_(parent_idx), child_idx_(child_idx), class_(INVALID_CLS), color_(c) {}

/* getters */
pt::id_t Edge::id() const { return this->id_; }
pt::idx_t Edge::get_class() const { return this->class_; }
pgt::color_e Edge::get_color() const { return this->color_; }
pt::idx_t Edge::get_parent_v_idx() const { return this->parent_idx_; }
pt::idx_t Edge::get_child_v_idx() const { return this->child_idx_; }
/* setters */
void Edge::set_class(pt::idx_t c) { this->class_ = c; }


/*
 * BackEdge
 * --------
 */
BackEdge::BackEdge(pt::id_t id,  pt::idx_t src,  pt::idx_t tgt, be_type_e t, pgt::color_e c)
  : id_(id), src(src), tgt(tgt), class_(INVALID_CLS), type_(t), color_(c) {}
/* getters */
pt::id_t BackEdge::id() const { return this->id_; }
pt::idx_t BackEdge::get_src() const { return this->src; }
pt::idx_t BackEdge::get_tgt() const { return this->tgt; }
pt::idx_t BackEdge::get_class() const { return this->class_; }
bool BackEdge::is_class_defined() const { return this->class_ != INVALID_CLS; }
be_type_e BackEdge::type() const { return this->type_; }
  //pgt::color_e BackEdge::get_color() const { return this->color_; }
/* setters */
void BackEdge::set_class(pt::idx_t c) { this->class_ = c; }

/*
 * Vertex
 * ------
 */

/* constructor(s) */
Vertex::Vertex(pt::idx_t dfs_num, pt::idx_t g_v_id, v_type_e type_)
    : dfs_num_(dfs_num), parent_e_idx_(pc::INVALID_IDX), hi_(pc::INVALID_IDX),
      g_v_id_(g_v_id), type_(type_) {}

// getters
pt::idx_t Vertex::g_v_id() const { return this->g_v_id_; }
v_type_e Vertex::type() const { return this->type_; }
pt::idx_t  Vertex::hi() const { return this->hi_; }
bool Vertex::is_root() const { return this->parent_e_idx_ == INVALID_IDX; }
bool Vertex::is_leaf() const { return this->child_e_idxs_.empty(); }
pt::idx_t Vertex::dfs_num() const { return this->dfs_num_; }
pt::idx_t Vertex::pre_order() const { return this->pre_order_; }
pt::idx_t Vertex::post_order() const { return this->post_order_; }

pt::idx_t Vertex::get_parent_e_idx() const { return this->parent_e_idx_; }

std::set<pt::idx_t> const &Vertex::get_ibe() const { return this->ibe; }
std::set<pt::idx_t> const &Vertex::get_obe() const { return this->obe; }

std::set<pt::idx_t> const &Vertex::get_child_edge_idxs() const { return this->child_e_idxs_; }
pt::idx_t Vertex::child_count() const { return static_cast<pt::idx_t>(this->child_e_idxs_.size()); }

// setters
void Vertex::add_obe(pt::idx_t obe_id) { this->obe.insert(obe_id); }
void Vertex::add_ibe(pt::idx_t ibe_id) { this->ibe.insert(ibe_id); }
void Vertex::add_child_e_idx(pt::idx_t e_id) { this->child_e_idxs_.insert(e_id); }
void Vertex::set_parent_e_idx(pt::idx_t e_idx) { this->parent_e_idx_ = e_idx; }
void Vertex::set_g_v_id(pt::idx_t g_v_id) { this->g_v_id_ = g_v_id; }
void Vertex::set_type(v_type_e t) { this->type_ = t; }
void Vertex::set_hi(pt::idx_t val) { this->hi_ = val; }
void Vertex::set_dfs_num(pt::idx_t idx) { this->dfs_num_ = idx; }
void Vertex::set_pre_order(pt::idx_t idx) { this->pre_order_ = idx; }
void Vertex::set_post_order(pt::idx_t idx) { this->post_order_ = idx; }

/*
 * Tree
 * ----
 */

// Constructor(s)

Tree::Tree(std::size_t size) :
  nodes(std::vector<Vertex>{}),
  tree_edges(std::vector<Edge>{}),
  back_edges(std::vector<BackEdge>{}),
  //bracket_lists(std::vector<BracketList*>{}),
  //sort_(std::vector<std::size_t>{}),
  //sort_g(std::vector<std::size_t>{}),
  equiv_class_count_(0) {
  this->nodes.reserve(size);
  this->tree_edges.reserve(size);
  this->back_edges.reserve(size);
  // this->bracket_lists.reserve(size);
  // TODO: is this necessary?
  this->bracket_lists = std::vector<WBracketList*>(size, nullptr);
}

Tree::~Tree(){
  this->nodes.clear();
  this->tree_edges.clear();
  this->back_edges.clear();
  for (std::size_t i = 0; i < this->bracket_lists.size(); ++i) {
    if (this->bracket_lists.at(i) != nullptr) {
      delete this->bracket_lists[i];
    }
    //else {
    //  std::cout << "Bracket list " << i << " is null" << std::endl;
    //}
  }
  this->bracket_lists.clear();

  //this->sort_.clear();
  //this->sort_g.clear();
}

  //void Tree::set_sort(std::size_t idx, std::size_t vertex) {
  //this->sort_.at(idx) = vertex;
  //}

  //void Tree::set_sort_g(std::size_t idx, std::size_t vertex) {
  //this->sort_g.at(idx) = vertex;
  //}

void Tree::set_dfs_num(std::size_t vertex, std::size_t dfs_num) {
  this->nodes.at(vertex).set_dfs_num(dfs_num);
}

void Tree::set_vertex_type(std::size_t vertex, v_type_e type) {
  this->nodes.at(vertex).set_type(type);
}

void Tree::add_vertex(Vertex&& v) {
  this->nodes.push_back(v);
}

Vertex& Tree::get_root()  { return this->nodes.at(this->get_root_idx()); }
std::size_t Tree::get_root_idx() const { return this->root_node_index; }

pt::idx_t Tree::vtx_count() const { return static_cast<pt::idx_t>(this->nodes.size()); };
  //pt::idx_t Tree::size() const { return this->nodes.size(); }
pt::idx_t Tree::tree_edge_count() const { return this->tree_edges.size(); }
pt::idx_t Tree::back_edge_count() const { return this->back_edges.size(); }

Vertex const &Tree::get_vertex(std::size_t vertex) const {
  return this->nodes.at(vertex);
}

Vertex& Tree::get_vertex_mut(std::size_t vertex) {
  return this->nodes.at(vertex);
}


Vertex const& Tree::get_p_vtx(std::size_t v_idx) const {
  std::size_t p_idx = this->get_parent_v_idx(v_idx);
  return this->get_vertex(p_idx);
}


std::size_t Tree::list_size(std::size_t vertex) {
  return this->bracket_lists.at(vertex)->size();
}

std::size_t Tree::get_hi(std::size_t vertex) {
  return this->nodes.at(vertex).hi();
}

bool Tree::is_desc(pt::idx_t a, pt::idx_t d) const {
  return this->get_vertex(a).pre_order() < this->get_vertex(d).pre_order() &&
         this->get_vertex(a).post_order() > this->get_vertex(d).post_order();
}

std::set<std::pair<std::size_t, std::size_t>>
Tree::get_children_w_id(std::size_t vertex) {
  std::set<std::pair<std::size_t, std::size_t>> res{};
  for (auto e_idx : this->nodes.at(vertex).get_child_edge_idxs()) {
    res.insert(
      std::make_pair(
        this->tree_edges.at(e_idx).id(),
        this->tree_edges.at(e_idx).get_child_v_idx())
   );
  }
  return res;
}

std::vector<Edge> Tree::get_child_edges(pt::idx_t v_idx) {
  return this->get_child_edges_mut(v_idx);
}

std::vector<Edge> Tree::get_child_edges_mut(pt::idx_t v_idx) {
  std::vector<Edge> v{};

  for (auto e_idx : this->nodes.at(v_idx).get_child_edge_idxs()) {
    v.push_back(this->tree_edges.at(e_idx));
  }

  return v;
}

std::vector<pt::idx_t> Tree::get_child_edge_idxs(pt::idx_t v_idx) const {
  std::vector<pt::idx_t> v {};
  for (auto e_idx : this->nodes.at(v_idx).get_child_edge_idxs()) {
    v.push_back(e_idx);
  }

  return v;
}

pt::idx_t Tree::get_child_count(pt::idx_t v_idx) const {
  return  this->nodes.at(v_idx).child_count();
}

Edge const& Tree::get_parent_edge(std::size_t vertex) const {
  return this->tree_edges.at(this->nodes.at(vertex).get_parent_e_idx());
}

std::set<std::size_t> Tree::get_obe_idxs(std::size_t vertex) const {
  std::set<std::size_t> v {};

  for (auto be_idx : this->nodes.at(vertex).get_obe()) {
    v.insert(be_idx);
  }

  return v;
}

std::set<std::size_t> Tree::get_ibe_idxs(std::size_t vertex) const {
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

std::set<std::pair<std::size_t, std::size_t>>
Tree::get_ibe_w_id(std::size_t vertex) {
  std::set<std::pair<std::size_t, std::size_t>> res{};
  for (auto e_idx : this->nodes.at(vertex).get_ibe()) {
    res.insert(
      std::make_pair(
        this->back_edges.at(e_idx).id(),
        this->back_edges.at(e_idx).get_src())
    );
  }
  return res;
}

std::set<std::size_t> Tree::get_children(pt::idx_t v_idx) const {
  std::set<std::size_t> res {};
  for (auto e_idx : this->nodes.at(v_idx).get_child_edge_idxs()) {
    res.insert(this->tree_edges.at(e_idx).get_child_v_idx());
  }
  return res;
}

std::set<std::size_t> Tree::get_ibe(std::size_t vertex) {
  std::set<std::size_t> res {};
  for (auto e_idx : this->nodes.at(vertex).get_ibe()) {
    res.insert(this->back_edges.at(e_idx).get_src());
  }

  return res;
}

std::set<pt::idx_t> Tree::get_ibe_src_v_idxs(std::size_t v_idx) const {
  std::set<pt::idx_t> res{};
  for (auto e_idx : this->nodes.at(v_idx).get_ibe()) {
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

std::set<pt::idx_t> Tree::get_obe_tgt_v_idxs(std::size_t v_idx) const {
  std::set<pt::idx_t> res{};
  for (auto e_idx : this->nodes.at(v_idx).get_obe()) {
    res.insert(this->back_edges.at(e_idx).get_tgt());
  }
  return res;
}

bool Tree::is_root(std::size_t vertex) const {
  return this->get_vertex(vertex).is_root();
}

bool Tree::is_leaf(std::size_t vertex) const {
  return this->get_vertex(vertex).is_leaf();
}

bool Tree::has_child(std::size_t vertex, std::size_t child_idx)  {
  return this->get_children(vertex).count(child_idx);
}

bool Tree::has_ibe(std::size_t vertex, std::size_t qry_idx)  {
  return this->get_ibe_src_v_idxs(vertex).count(qry_idx);
}

bool Tree::has_obe(std::size_t vertex, std::size_t qry_idx)  {
  return this->get_obe_tgt_v_idxs(vertex).count(qry_idx);
}

Edge& Tree::get_incoming_edge(std::size_t vertex) {
  std::size_t e_idx = this->nodes.at(vertex).get_parent_e_idx();
  return this->tree_edges.at(e_idx);
}

std::size_t Tree::get_parent(std::size_t vertex) {
  std::size_t e_idx = this->get_vertex(vertex).get_parent_e_idx();
  return this->tree_edges.at(e_idx).get_parent_v_idx();
}

std::size_t Tree::get_parent_v_idx(std::size_t v_idx) const {
  std::size_t e_idx = this->get_vertex(v_idx).get_parent_e_idx();
  return this->tree_edges.at(e_idx).get_parent_v_idx();
}

const Edge& Tree::get_tree_edge(std::size_t edge_idx) const {
    return this->tree_edges.at(edge_idx);
}

std::size_t Tree::get_graph_edge_id(std::size_t tree_edge_id) const {
    return this->tree_graph_idx_map_.at(tree_edge_id);
}

BackEdge& Tree::get_backedge(std::size_t backedge_idx) {
  return this->back_edges.at(backedge_idx);
}

const BackEdge &Tree::get_be(std::size_t backedge_idx) const {
  return this->back_edges.at(backedge_idx);
}

BackEdge &Tree::get_backedge_ref_given_id(std::size_t backedge_id) {
  std::size_t be_idx = this->be_id_to_idx_map_.at(backedge_id);
  return this->back_edges.at(be_idx);
}

BackEdge Tree::get_backedge_given_id(std::size_t backedge_id) {
  std::size_t be_idx = this->be_id_to_idx_map_.at(backedge_id);
  return this->back_edges[be_idx];
}

void Tree::add_tree_edge(pt::idx_t frm, pt::idx_t to, pgt::color_e c) {
  std::size_t edge_idx = this->tree_edges.size();
  std::size_t edge_count = edge_idx + this->back_edges.size();
  this->tree_edges.push_back(Edge(edge_count, frm, to, c));

  this->nodes[frm].add_child_e_idx(edge_idx);
  this->nodes[to].set_parent_e_idx(edge_idx);
}

pt::idx_t Tree::add_be(pt::idx_t frm, pt::idx_t to, be_type_e t, pgt::color_e c) {
  pt::idx_t back_edge_idx = this->back_edges.size();
  pt::idx_t edge_count = back_edge_idx + this->tree_edges.size();
  this->back_edges.push_back(BackEdge(edge_count, frm, to, t, c));
  this->nodes[frm].add_obe(back_edge_idx);
  this->nodes[to].add_ibe(back_edge_idx);

  this->be_id_to_idx_map_[edge_count] = back_edge_idx;

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
void Tree::concat_bracket_lists(std::size_t parent_vertex, std::size_t child_vertex) {
  std::string fn_name = std::format("[povu::spanning_tree::Tree::{}]", __func__);

  WBracketList* bl_p = this->bracket_lists[parent_vertex];
  WBracketList* bl_c = this->bracket_lists[child_vertex];

  if (bl_p == nullptr) {
    this->bracket_lists[parent_vertex] = bl_c;
    this->bracket_lists[child_vertex] = nullptr;
  }
  else {
    bl_p->concat(bl_c);
  }
}

// TODO: once deleted do we care to reflect changes in the concated ones?
// WE DO!!
/**
 * delete the bracket that is associated with the backedge
 * given a vertex id and a backedge idx
 */
void Tree::del_bracket(std::size_t vertex, std::size_t backedge_idx) {
  std::string fn_name = std::format("[povu::spanning_tree::Tree::{}]", __func__);

  std::size_t be_id = this->back_edges.at(backedge_idx).id();
  this->bracket_lists[vertex]->del(be_id);
}


void Tree::push(std::size_t vertex, std::size_t backege_idx) {
  std::string fn_name = std::format("[povu::spanning_tree::{}]", __func__);

  // TODO: based on the Tree constructor we expect the pointer at v_idx will
  // never be null why then do we need to check for null else code fails
  // we then create a  bracket using the backedge ID

  if (this->bracket_lists[vertex] == nullptr) {
    this->bracket_lists[vertex] = new WBracketList{};
  }

  this->bracket_lists[vertex]->push(Bracket(this->back_edges.at(backege_idx).id()));
}


BracketList& Tree::get_bracket_list(std::size_t vertex) {
  std::string fn_name = std::format("[povu::spanning_tree::{}]", __func__);
  if (this->bracket_lists[vertex] == nullptr) {
    throw std::runtime_error(std::format("{} Bracket list is null", fn_name));
  }

  return this->bracket_lists[vertex]->get_bracket_list();
}


Bracket& Tree::top(std::size_t vertex) {
  std::string fn_name = std::format("[povu::spanning_tree::{}]", __func__);

  return this->bracket_lists[vertex]->top();
}


std::size_t Tree::new_class() { return this->equiv_class_count_++; }


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
//std::size_t Tree::get_sorted(std::size_t idx) { return  this->sort_.at(idx);}

//std::size_t Tree::get_sorted_g(std::size_t idx) { return  this->sort_g.at(idx);}

const std::map<std::size_t, std::pair<be_type_e, std::size_t>>& Tree::get_g_edge_idx_map() const {
  return this->g_edge_idx_map;
}

void Tree::print_dot(std::ostream &os) {

  /* ---------- Helper Functions ---------- */

  auto vtx_to_dot = [&](std::size_t i) {
    const Vertex &vertex = this->get_vertex(i);
    std::string str;

    switch (vertex.type()) {
    case v_type_e::dummy:
      str = std::format("\t{} [style=filled, fillcolor=pink];\n", i);
      break;
    case v_type_e::l:
    case v_type_e::r:
      std::string sign = (vertex.type() == pgt::v_type_e::l) ? "+" : "-";
      str = std::format(
          "\t{} [style=filled, fillcolor=lightblue, label = \"{} \\n ({}{}) \\n [{},{}]\"];\n",
                                                                                            i, i, vertex.g_v_id(), sign, vertex.pre_order(), vertex.post_order());
      break;
    }

    os << str;
  };

  auto tree_edge_to_dot = [&](pt::idx_t p_v_idx, Edge &e) {
    std::string cls = e.get_class() == INVALID_CLS ? "" : std::to_string(e.get_class());
    std::string clr = e.get_color() == pgt::color_e::gray ? "gray" : "black";

    os << std::format("\t{}  -- {}  [label=\"{} {}\" color={}];\n",
                             p_v_idx, e.get_child_v_idx(), e.id(), cls, clr);
  };

  auto be_to_dot = [&](pt::idx_t i, const std::pair<std::size_t, std::size_t> &o) {
    auto [f,s]= o;
    std::string cl = f > 10000 ? "\u2205" : std::to_string(f);
    // a capping backedge is red and can never have been gray
    BackEdge be = this->get_backedge_given_id(f);
    std::string class_ = be.get_class() == INVALID_CLS ? "" : std::to_string(be.get_class());

    std::string color = [&]() -> std::string {
      switch (be.type()) {
      case be_type_e::capping_back_edge: return pc::RED;
      case be_type_e::simplifying_back_edge: return pc::BLUE;
      default: return pc::GRAY; // "Normal" backedge
      }
    }();

    os << std::format("\t{} -- {} [label=\"{} {}\" style=\"dotted\" "
                             "penwidth=\"3\" color=\"{}\"];\n",
                             i, be.get_tgt(), cl, class_, color);
  };

  /* ---------- dot format header ---------- */

  os << std::format(
    "graph G {{\n"
    "\trankdir = LR;\n"
    "\tnode[shape = circle];\n"
    "\tedge [arrowhead=vee];\n"
  );

  //print the vertices
  for (std::size_t i{}; i < this->vtx_count(); i++){
    vtx_to_dot(i);
  }

  // print the edges
  for (std::size_t i{}; i < this->vtx_count(); i++) {
    for (auto &c : this->get_child_edges_mut(i)) { // tree edges
      tree_edge_to_dot(i, c);
    }

    for (auto o : this->get_obe_w_id(i)) { // back edges
      be_to_dot(i, o);
    }
  }

  // end the dot format
  os << "}" << std::endl;
}

} // namespace spanning_tree
