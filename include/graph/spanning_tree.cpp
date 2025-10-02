#include "./spanning_tree.hpp"

namespace povu::spanning_tree
{

/*
 * Edge
 * ----
 */
/*  constructor(s) */
Edge::Edge(pt::id_t id, pt::idx_t parent_idx, pt::idx_t child_idx,
	   pgt::color_e c)
    : id_(id), parent_idx_(parent_idx), child_idx_(child_idx),
      class_(INVALID_CLS), color_(c)
{}

/* getters */
pt::id_t Edge::id() const
{
	return this->id_;
}

pt::idx_t Edge::get_class() const
{
	return this->class_;
}

pgt::color_e Edge::get_color() const
{
	return this->color_;
}

pt::idx_t Edge::get_parent_v_idx() const
{
	return this->parent_idx_;
}

pt::idx_t Edge::get_child_v_idx() const
{
	return this->child_idx_;
}

/* setters */
void Edge::set_class(pt::idx_t c)
{
	this->class_ = c;
}

/*
 * BackEdge
 * --------
 */
BackEdge::BackEdge(pt::id_t id, pt::idx_t src, pt::idx_t tgt, be_type_e t,
		   pgt::color_e c)
    : id_(id), src_(src), tgt_(tgt), class_(INVALID_CLS), type_(t), color_(c)
{}

/* getters */
pt::id_t BackEdge::id() const
{
	return this->id_;
}

pt::idx_t BackEdge::get_src() const
{
	return this->src_;
}

pt::idx_t BackEdge::get_tgt() const
{
	return this->tgt_;
}

pt::idx_t BackEdge::get_class() const
{
	return this->class_;
}

bool BackEdge::is_class_defined() const
{
	return this->class_ != INVALID_CLS;
}

be_type_e BackEdge::type() const
{
	return this->type_;
}

/* setters */
void BackEdge::set_class(pt::idx_t c)
{
	this->class_ = c;
}

// ======
// Vertex
// ======

/* constructor(s) */
Vertex::Vertex(pt::idx_t dfs_num, pt::idx_t g_v_id, v_type_e type)
    : dfs_num_(dfs_num), parent_e_idx_(pc::INVALID_IDX), hi_(pc::INVALID_IDX),
      g_v_id_(g_v_id), pre_order_(pc::INVALID_IDX),
      post_order_(pc::INVALID_IDX), type_(type)
{}

// getters
pt::idx_t Vertex::g_v_id() const
{
	return this->g_v_id_;
}

v_type_e Vertex::type() const
{
	return this->type_;
}

pt::idx_t Vertex::hi() const
{
	return this->hi_;
}

bool Vertex::is_root() const
{
	return this->parent_e_idx_ == INVALID_IDX;
}

bool Vertex::is_leaf() const
{
	return this->child_e_idxs_.empty();
}

pt::idx_t Vertex::dfs_num() const
{
	return this->dfs_num_;
}

pt::idx_t Vertex::pre_order() const
{
	return this->pre_order_;
}

pt::idx_t Vertex::post_order() const
{
	return this->post_order_;
}

pt::idx_t Vertex::get_parent_e_idx() const
{
	return this->parent_e_idx_;
}

std::set<pt::idx_t> const &Vertex::get_ibe() const
{
	return this->ibe;
}

std::set<pt::idx_t> const &Vertex::get_obe() const
{
	return this->obe;
}

std::set<pt::idx_t> const &Vertex::get_child_edge_idxs() const
{
	return this->child_e_idxs_;
}

pt::idx_t Vertex::child_count() const
{
	return static_cast<pt::idx_t>(this->child_e_idxs_.size());
}

// setters
void Vertex::add_obe(pt::idx_t obe_id)
{
	this->obe.insert(obe_id);
}

void Vertex::add_ibe(pt::idx_t ibe_id)
{
	this->ibe.insert(ibe_id);
}

void Vertex::add_child_e_idx(pt::idx_t e_id)
{
	this->child_e_idxs_.insert(e_id);
}

void Vertex::set_parent_e_idx(pt::idx_t e_idx)
{
	this->parent_e_idx_ = e_idx;
}

void Vertex::set_g_v_id(pt::idx_t g_v_id)
{
	this->g_v_id_ = g_v_id;
}

void Vertex::set_type(v_type_e t)
{
	this->type_ = t;
}

void Vertex::set_hi(pt::idx_t val)
{
	this->hi_ = val;
}

void Vertex::set_dfs_num(pt::idx_t idx)
{
	this->dfs_num_ = idx;
}

void Vertex::set_pre_order(pt::idx_t idx)
{
	this->pre_order_ = idx;
}

void Vertex::set_post_order(pt::idx_t idx)
{
	this->post_order_ = idx;
}

/*
 * Tree
 * ----
 */

// Constructor(s)

Tree::Tree(std::size_t size)
    : nodes(std::vector<Vertex>{}), tree_edges(std::vector<Edge>{}),
      back_edges(std::vector<BackEdge>{}),
      // bracket_lists(std::vector<BracketList*>{}),
      // sort_(std::vector<std::size_t>{}),
      // sort_g(std::vector<std::size_t>{}),
      equiv_class_count_(0)
{
	this->nodes.reserve(size);
	this->tree_edges.reserve(size);
	this->back_edges.reserve(size);
	// this->bracket_lists.reserve(size);
	// TODO: is this necessary?
	this->bracket_lists = std::vector<WBracketList *>(size, nullptr);
}

Tree Tree::from_bd(const bd::VG &g)
{

	const bool has_tips{!g.tips().empty()};
	const pt::idx_t root_idx{0};

	std::stack<pt::idx_t> s;
	std::vector<u_int8_t> visited(g.vtx_count(), 0);

	std::set<std::pair<pt::idx_t, pt::idx_t>> added_edges;
	std::unordered_set<pt::idx_t> self_loops;

	pt::idx_t order{};
	// pt::idx_t post_order {};

	pt::idx_t counter{0};		 // dfs num
	bool found_new_neighbour{false}; // neighbours exhausted

	pt::idx_t p_idx{pc::INVALID_IDX}; // parent_idx

	pt::idx_t t_vtx_count =
		has_tips ? (2 * g.vtx_count()) + 1 : 2 * g.vtx_count();
	Tree t{t_vtx_count};

	// biedged idx to tree idx (or counter)
	std::vector<pt::id_t> be_idx_to_ctr(t_vtx_count, 0);

	auto unordered_pair = [](pt::idx_t a,
				 pt::idx_t b) -> std::pair<pt::idx_t, pt::idx_t>
	{
		return {std::min(a, b), std::max(a, b)};
	};

	auto connect = [&](pt::idx_t a, pt::idx_t b) -> void
	{
		added_edges.insert(unordered_pair(a, b));
	};

	auto are_connected = [&](pt::idx_t a, pt::idx_t b) -> bool
	{
		return pv_cmp::contains(added_edges, unordered_pair(a, b));
	};

	auto to_be = [&g](pgt::side_n_id_t i) -> pt::idx_t
	{
		auto [ve, v_id] = i;
		pt::idx_t v_idx = g.v_id_to_idx(v_id);
		return (ve == pgt::v_end_e::l) ? v_idx * 2 : (v_idx * 2) + 1;
	};

	auto to_bd = [&g](pt::idx_t be_v_idx) -> pgt::side_n_idx_t
	{
		pgt::v_end_e ve =
			(be_v_idx % 2 == 0) ? pgt::v_end_e::l : pgt::v_end_e::r;
		pt::id_t v_id = g.v_idx_to_id(be_v_idx / 2);
		return {ve, v_id};
	};

	auto end2typ = [](pgt::v_end_e e) -> pgt::v_type_e
	{
		return pgt::v_end_e::l == e ? pgt::v_type_e::l
					    : pgt::v_type_e::r;
	};

	auto add_vertex_to_tree = [&](pgt::v_end_e e,
				      pt::idx_t bd_v_idx) -> void
	{
		const bd::Vertex &v = g.get_vertex_by_idx(bd_v_idx);

		Vertex v1{counter++, v.id(), end2typ(e)};
		v1.set_pre_order(order++);
		t.add_vertex(std::move(v1));

		Vertex v2{counter++, v.id(), end2typ(pgt::complement(e))};
		v2.set_pre_order(order++);
		t.add_vertex(std::move(v2));
		// t.add_vertex({counter++, v.id(),
		// end2typ(pgt::complement(e))});

		be_idx_to_ctr[to_be({e, v.id()})] = counter - 2;
		be_idx_to_ctr[to_be({pgt::complement(e), v.id()})] =
			counter - 1;

		// add edges
		if (p_idx != pc::INVALID_IDX) {
			t.add_tree_edge(p_idx, counter - 2, pgt::color_e::gray);
			connect(p_idx, counter - 2);
		}
		t.add_tree_edge(counter - 2, counter - 1, pgt::color_e::black);
		connect(counter - 2, counter - 1);
	};

	// returns true if it discovers a new vertex (neighbour), false
	// otherwise
	auto process_edge = [&](pt::idx_t bd_v_idx, pgt::v_end_e ve,
				pt::idx_t e_idx) -> bool
	{
		const bd::Edge &e = g.get_edge(e_idx);

		auto [os, ov_idx] =
			e.get_other_vtx(bd_v_idx, ve); // o for other
		pt::idx_t o_be_idx = to_be({os, g.v_idx_to_id(ov_idx)});

		if (!visited[ov_idx]) { // has not been visited
			add_vertex_to_tree(os, ov_idx);

			visited[ov_idx] = 1;
			s.push(to_be({os, g.v_idx_to_id(ov_idx)}));
			s.push(to_be(
				{pgt::complement(os), g.v_idx_to_id(ov_idx)}));

			return true;
		}
		else if (!are_connected(p_idx, be_idx_to_ctr[o_be_idx])) {
			// add a backedge if:
			//  - not a parent child relationship
			//  - a backedge does not already exist
			t.add_be(p_idx, be_idx_to_ctr[o_be_idx],
				 be_type_e::back_edge, pgt::color_e::gray);
			connect(p_idx, be_idx_to_ctr[o_be_idx]);
		}
		else if (__builtin_expect(
				 (bd_v_idx == ov_idx &&
				  !pv_cmp::contains(self_loops, bd_v_idx)),
				 0)) {
			// add a self loop backedge, a parent-child relationship
			t.add_be(p_idx, be_idx_to_ctr[o_be_idx],
				 be_type_e::back_edge, pgt::color_e::gray);
			self_loops.insert(bd_v_idx);
		}

		return false;
	};

	if (has_tips) { // add a dummy vertex to the tree
		p_idx = counter;
		Vertex v{counter++, pc::DUMMY_VTX_ID, v_type_e::dummy};
		v.set_pre_order(order++);
		t.add_vertex(std::move(v));
	}

	side_n_id_t start =
		has_tips ? *g.tips().begin()
			 : pgt::side_n_id_t{pgt::v_end_e::l, g.v_idx_to_id(0)};
	auto [s_v_end, s_v_id] = start;
	pt::idx_t s_v_idx = g.v_id_to_idx(s_v_id);
	s.push(to_be({s_v_end, s_v_id}));
	s.push(to_be({pgt::complement(s_v_end), s_v_id}));
	visited[s_v_idx] = 1;
	add_vertex_to_tree(s_v_end, s_v_idx);

	/* ---------- Main Loop ---------- */

	while (!s.empty()) {
		found_new_neighbour = false;
		pt::idx_t be_v_idx = s.top();

		p_idx = be_idx_to_ctr[be_v_idx];
		auto [syd, v_id] = to_bd(be_v_idx);
		pt::idx_t bd_v_idx = g.v_id_to_idx(v_id);

		const bd::Vertex &v = g.get_vertex_by_id(v_id);
		const std::set<pt::idx_t> &neighbours =
			syd == pgt::v_end_e::l ? v.get_edges_l()
					       : v.get_edges_r();

		// if no neighbours then it is a tip. Add a backedge to the root
		if (__builtin_expect((neighbours.empty() &&
				      !are_connected(p_idx, root_idx)),
				     0)) {
			t.add_be(p_idx, root_idx, be_type_e::back_edge,
				 pgt::color_e::gray);
			connect(p_idx, root_idx);
		}

		// stop processing edges on the first instance of finding a new
		// neighbour
		for (auto e_idx : neighbours) {
			if ((found_new_neighbour =
				     process_edge(bd_v_idx, syd, e_idx))) {
				break;
			}
		}

		if (!found_new_neighbour) {
			t.get_vertex_mut(p_idx).set_post_order(order++);
			s.pop();
		}
	}

	if (has_tips) {
		t.get_vertex_mut(0).set_post_order(order++);
	}

	return t;
}

Tree::~Tree()
{
	this->nodes.clear();
	this->tree_edges.clear();
	this->back_edges.clear();
	for (std::size_t i = 0; i < this->bracket_lists.size(); ++i) {
		if (this->bracket_lists.at(i) != nullptr) {
			delete this->bracket_lists[i];
		}
		// else {
		//   std::cout << "Bracket list " << i << " is null" <<
		//   std::endl;
		// }
	}
	this->bracket_lists.clear();

	// this->sort_.clear();
	// this->sort_g.clear();
}

// void Tree::set_sort(std::size_t idx, std::size_t vertex) {
// this->sort_.at(idx) = vertex;
// }

// void Tree::set_sort_g(std::size_t idx, std::size_t vertex) {
// this->sort_g.at(idx) = vertex;
// }

void Tree::set_dfs_num(std::size_t vertex, std::size_t dfs_num)
{
	this->nodes.at(vertex).set_dfs_num(dfs_num);
}

void Tree::set_vertex_type(std::size_t vertex, v_type_e type)
{
	this->nodes.at(vertex).set_type(type);
}

void Tree::add_vertex(Vertex &&v)
{
	this->nodes.push_back(v);
}

Vertex &Tree::get_root()
{
	return this->nodes.at(this->get_root_idx());
}

std::size_t Tree::get_root_idx() const
{
	return this->root_node_index;
}

pt::idx_t Tree::vtx_count() const
{
	return static_cast<pt::idx_t>(this->nodes.size());
};

// pt::idx_t Tree::size() const { return this->nodes.size(); }
pt::idx_t Tree::tree_edge_count() const
{
	return this->tree_edges.size();
}

pt::idx_t Tree::back_edge_count() const
{
	return this->back_edges.size();
}

Vertex const &Tree::get_vertex(std::size_t vertex) const
{
	return this->nodes.at(vertex);
}

Vertex &Tree::get_vertex_mut(std::size_t vertex)
{
	return this->nodes.at(vertex);
}

Vertex const &Tree::get_p_vtx(std::size_t v_idx) const
{
	std::size_t p_idx = this->get_parent_v_idx(v_idx);
	return this->get_vertex(p_idx);
}

std::size_t Tree::list_size(std::size_t vertex)
{
	return this->bracket_lists.at(vertex)->size();
}

std::size_t Tree::get_hi(std::size_t vertex)
{
	return this->nodes.at(vertex).hi();
}

bool Tree::is_desc(pt::idx_t a, pt::idx_t d) const
{
	return this->get_vertex(a).pre_order() <
		       this->get_vertex(d).pre_order() &&
	       this->get_vertex(a).post_order() >
		       this->get_vertex(d).post_order();
}

std::set<std::pair<std::size_t, std::size_t>>
Tree::get_children_w_id(std::size_t vertex)
{
	std::set<std::pair<std::size_t, std::size_t>> res{};
	for (auto e_idx : this->nodes.at(vertex).get_child_edge_idxs()) {
		res.insert(std::make_pair(
			this->tree_edges.at(e_idx).id(),
			this->tree_edges.at(e_idx).get_child_v_idx()));
	}
	return res;
}

std::vector<Edge> Tree::get_child_edges(pt::idx_t v_idx)
{
	return this->get_child_edges_mut(v_idx);
}

std::vector<Edge> Tree::get_child_edges_mut(pt::idx_t v_idx)
{
	std::vector<Edge> v{};

	for (auto e_idx : this->nodes.at(v_idx).get_child_edge_idxs()) {
		v.push_back(this->tree_edges.at(e_idx));
	}

	return v;
}

std::vector<pt::idx_t> Tree::get_child_edge_idxs(pt::idx_t v_idx) const
{
	std::vector<pt::idx_t> v{};
	for (auto e_idx : this->nodes.at(v_idx).get_child_edge_idxs()) {
		v.push_back(e_idx);
	}

	return v;
}

pt::idx_t Tree::get_child_count(pt::idx_t v_idx) const
{
	return this->nodes.at(v_idx).child_count();
}

Edge const &Tree::get_parent_edge(std::size_t vertex) const
{
	return this->tree_edges.at(this->nodes.at(vertex).get_parent_e_idx());
}

std::set<std::size_t> Tree::get_obe_idxs(std::size_t vertex) const
{
	std::set<std::size_t> v{};

	for (auto be_idx : this->nodes.at(vertex).get_obe()) {
		v.insert(be_idx);
	}

	return v;
}

std::set<std::size_t> Tree::get_ibe_idxs(std::size_t vertex) const
{
	std::set<std::size_t> v{};

	for (auto be_idx : this->nodes.at(vertex).get_ibe()) {
		v.insert(be_idx);
	}

	return v;
}

std::set<std::pair<std::size_t, std::size_t>>
Tree::get_obe_w_id(std::size_t vertex)
{
	std::set<std::pair<std::size_t, std::size_t>> res{};
	for (auto e_idx : this->nodes.at(vertex).get_obe()) {
		res.insert(
			std::make_pair(this->back_edges.at(e_idx).id(),
				       this->back_edges.at(e_idx).get_tgt()));
	}
	return res;
}

std::set<std::pair<std::size_t, std::size_t>>
Tree::get_ibe_w_id(std::size_t vertex)
{
	std::set<std::pair<std::size_t, std::size_t>> res{};
	for (auto e_idx : this->nodes.at(vertex).get_ibe()) {
		res.insert(
			std::make_pair(this->back_edges.at(e_idx).id(),
				       this->back_edges.at(e_idx).get_src()));
	}
	return res;
}

std::set<std::size_t> Tree::get_children(pt::idx_t v_idx) const
{
	std::set<std::size_t> res{};
	for (auto e_idx : this->nodes.at(v_idx).get_child_edge_idxs()) {
		res.insert(this->tree_edges.at(e_idx).get_child_v_idx());
	}
	return res;
}

std::set<std::size_t> Tree::get_ibe(std::size_t vertex)
{
	std::set<std::size_t> res{};
	for (auto e_idx : this->nodes.at(vertex).get_ibe()) {
		res.insert(this->back_edges.at(e_idx).get_src());
	}

	return res;
}

std::set<pt::idx_t> Tree::get_ibe_src_v_idxs(std::size_t v_idx) const
{
	std::set<pt::idx_t> res{};
	for (auto e_idx : this->nodes.at(v_idx).get_ibe()) {
		res.insert(this->back_edges.at(e_idx).get_src());
	}

	return res;
}

std::set<std::size_t> Tree::get_obe(std::size_t vertex)
{
	std::set<std::size_t> res{};
	for (auto e_idx : this->nodes.at(vertex).get_obe()) {
		res.insert(this->back_edges.at(e_idx).get_tgt());
	}
	return res;
}

std::set<pt::idx_t> Tree::get_obe_tgt_v_idxs(std::size_t v_idx) const
{
	std::set<pt::idx_t> res{};
	for (auto e_idx : this->nodes.at(v_idx).get_obe()) {
		res.insert(this->back_edges.at(e_idx).get_tgt());
	}
	return res;
}

bool Tree::is_root(std::size_t vertex) const
{
	return this->get_vertex(vertex).is_root();
}

bool Tree::is_leaf(std::size_t vertex) const
{
	return this->get_vertex(vertex).is_leaf();
}

bool Tree::has_child(std::size_t vertex, std::size_t child_idx)
{
	return this->get_children(vertex).count(child_idx);
}

bool Tree::has_ibe(std::size_t vertex, std::size_t qry_idx)
{
	return this->get_ibe_src_v_idxs(vertex).count(qry_idx);
}

bool Tree::has_obe(std::size_t vertex, std::size_t qry_idx)
{
	return this->get_obe_tgt_v_idxs(vertex).count(qry_idx);
}

Edge &Tree::get_incoming_edge(std::size_t vertex)
{
	std::size_t e_idx = this->nodes.at(vertex).get_parent_e_idx();
	return this->tree_edges.at(e_idx);
}

std::size_t Tree::get_parent(std::size_t vertex)
{
	std::size_t e_idx = this->get_vertex(vertex).get_parent_e_idx();
	return this->tree_edges.at(e_idx).get_parent_v_idx();
}

std::size_t Tree::get_parent_v_idx(std::size_t v_idx) const
{
	std::size_t e_idx = this->get_vertex(v_idx).get_parent_e_idx();
	return this->tree_edges.at(e_idx).get_parent_v_idx();
}

const Edge &Tree::get_tree_edge(std::size_t edge_idx) const
{
	return this->tree_edges.at(edge_idx);
}

std::size_t Tree::get_graph_edge_id(std::size_t tree_edge_id) const
{
	return this->tree_graph_idx_map_.at(tree_edge_id);
}

BackEdge &Tree::get_backedge(std::size_t backedge_idx)
{
	return this->back_edges.at(backedge_idx);
}

const BackEdge &Tree::get_be(std::size_t backedge_idx) const
{
	return this->back_edges.at(backedge_idx);
}

BackEdge &Tree::get_backedge_ref_given_id(std::size_t backedge_id)
{
	std::size_t be_idx = this->be_id_to_idx_map_.at(backedge_id);
	return this->back_edges.at(be_idx);
}

BackEdge Tree::get_backedge_given_id(std::size_t backedge_id)
{
	std::size_t be_idx = this->be_id_to_idx_map_.at(backedge_id);
	return this->back_edges[be_idx];
}

void Tree::add_tree_edge(pt::idx_t frm, pt::idx_t to, pgt::color_e c)
{
	std::size_t edge_idx = this->tree_edges.size();
	std::size_t edge_count = edge_idx + this->back_edges.size();
	this->tree_edges.push_back(Edge(edge_count, frm, to, c));

	this->nodes[frm].add_child_e_idx(edge_idx);
	this->nodes[to].set_parent_e_idx(edge_idx);
}

pt::idx_t Tree::add_be(pt::idx_t frm, pt::idx_t to, be_type_e t, pgt::color_e c)
{
	pt::idx_t back_edge_idx = this->back_edges.size();
	pt::idx_t edge_count = back_edge_idx + this->tree_edges.size();
	this->back_edges.push_back(BackEdge(edge_count, frm, to, t, c));
	this->nodes[frm].add_obe(back_edge_idx);
	this->nodes[to].add_ibe(back_edge_idx);

	this->be_id_to_idx_map_[edge_count] = back_edge_idx;

	return back_edge_idx;
}

void Tree::set_hi(std::size_t vertex, std::size_t val)
{
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
void Tree::concat_bracket_lists(std::size_t parent_vertex,
				std::size_t child_vertex)
{
	std::string fn_name =
		pv_cmp::format("[povu::spanning_tree::Tree::{}]", __func__);

	WBracketList *bl_p = this->bracket_lists[parent_vertex];
	WBracketList *bl_c = this->bracket_lists[child_vertex];

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
void Tree::del_bracket(std::size_t vertex, std::size_t backedge_idx)
{
	std::string fn_name =
		pv_cmp::format("[povu::spanning_tree::Tree::{}]", __func__);

	std::size_t be_id = this->back_edges.at(backedge_idx).id();
	this->bracket_lists[vertex]->del(be_id);
}

void Tree::push(std::size_t vertex, std::size_t backege_idx)
{
	std::string fn_name =
		pv_cmp::format("[povu::spanning_tree::{}]", __func__);

	// TODO: based on the Tree constructor we expect the pointer at v_idx
	// will never be null why then do we need to check for null else code
	// fails we then create a  bracket using the backedge ID

	if (this->bracket_lists[vertex] == nullptr) {
		this->bracket_lists[vertex] = new WBracketList{};
	}

	this->bracket_lists[vertex]->push(
		Bracket(this->back_edges.at(backege_idx).id()));
}

BracketList &Tree::get_bracket_list(std::size_t vertex)
{
	std::string fn_name =
		pv_cmp::format("[povu::spanning_tree::{}]", __func__);
	if (this->bracket_lists[vertex] == nullptr) {
		throw std::runtime_error(
			pv_cmp::format("{} Bracket list is null", fn_name));
	}

	return this->bracket_lists[vertex]->get_bracket_list();
}

Bracket &Tree::top(std::size_t vertex)
{
	std::string fn_name =
		pv_cmp::format("[povu::spanning_tree::{}]", __func__);

	return this->bracket_lists[vertex]->top();
}

std::size_t Tree::new_class()
{
	return this->equiv_class_count_++;
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
// std::size_t Tree::get_sorted(std::size_t idx) { return  this->sort_.at(idx);}

// std::size_t Tree::get_sorted_g(std::size_t idx) { return
// this->sort_g.at(idx);}

const std::map<std::size_t, std::pair<be_type_e, std::size_t>> &
Tree::get_g_edge_idx_map() const
{
	return this->g_edge_idx_map;
}

void Tree::print_dot(std::ostream &os)
{

	/* ---------- Helper Functions ---------- */

	auto vtx_to_dot = [&](std::size_t i)
	{
		const Vertex &vertex = this->get_vertex(i);
		std::string str;

		switch (vertex.type()) {
		case v_type_e::dummy:
			str = pv_cmp::format(
				"\t{} [style=filled, fillcolor=pink];\n", i);
			break;
		case v_type_e::l:
		case v_type_e::r:
			std::string sign =
				(vertex.type() == pgt::v_type_e::l) ? "+" : "-";
			str = pv_cmp::format(
				"\t{} [style=filled, fillcolor=lightblue, "
				"label = \"{} \\n ({}{}) \\n [{},{}]\"];\n",
				i, i, vertex.g_v_id(), sign, vertex.pre_order(),
				vertex.post_order());
			break;
		}

		os << str;
	};

	auto tree_edge_to_dot = [&](pt::idx_t p_v_idx, Edge &e)
	{
		std::string cls = e.get_class() == INVALID_CLS
					  ? ""
					  : std::to_string(e.get_class());
		std::string clr =
			e.get_color() == pgt::color_e::gray ? "gray" : "black";

		os << pv_cmp::format(
			"\t{}  -- {}  [label=\"{} {}\" color={}];\n", p_v_idx,
			e.get_child_v_idx(), e.id(), cls, clr);
	};

	auto be_to_dot =
		[&](pt::idx_t i, const std::pair<std::size_t, std::size_t> &o)
	{
		auto [f, s] = o;
		std::string cl = f > 10000 ? "\u2205" : std::to_string(f);
		// a capping backedge is red and can never have been gray
		BackEdge be = this->get_backedge_given_id(f);
		std::string class_ = be.get_class() == INVALID_CLS
					     ? ""
					     : std::to_string(be.get_class());

		std::string_view color = [&]() -> std::string_view
		{
			switch (be.type()) {
			case be_type_e::capping_back_edge:
				return pc::RED;
			case be_type_e::simplifying_back_edge:
				return pc::BLUE;
			default:
				return pc::GRAY; // "Normal" backedge
			}
		}();

		os << pv_cmp::format(
			"\t{} -- {} [label=\"{} {}\" style=\"dotted\" "
			"penwidth=\"3\" color=\"{}\"];\n",
			i, be.get_tgt(), cl, class_, color);
	};

	/* ---------- dot format header ---------- */

	os << "graph G {\n"
	      "\trankdir = LR;\n"
	      "\tnode[shape = circle];\n"
	      "\tedge [arrowhead=vee];\n";

	// print the vertices
	for (std::size_t i{}; i < this->vtx_count(); i++)
		vtx_to_dot(i);

	// print the edges
	for (std::size_t i{}; i < this->vtx_count(); i++) {
		for (auto &c : this->get_child_edges_mut(i)) // tree edges
			tree_edge_to_dot(i, c);

		for (auto o : this->get_obe_w_id(i)) // back edges
			be_to_dot(i, o);
	}

	// end the dot format
	os << "}" << std::endl;
}

} // namespace povu::spanning_tree
