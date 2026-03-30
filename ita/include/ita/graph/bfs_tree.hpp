#ifndef IT_BFS_TREE_HPP
#define IT_BFS_TREE_HPP

#include <list>
#include <ostream> // for ostream
#include <queue>   // for queue
#include <string>
#include <vector> // for vector

#include "povu/common/compat.hpp"    // for pv_cmp, format
#include "povu/common/constants.hpp" // for INVALID_IDX
#include "povu/common/core.hpp"	     // for pt
#include "povu/graph/types.hpp"	     // for id_or_t, op_t

namespace ita::bfs
{

inline constexpr std::string_view MODULE = "povu::graph::bfs";

enum edge_type : pt::u8 {
	tree_edge,
	cross_edge
};

// Directed edge
struct directed_edge {
	pt::u32 from;
	pt::u32 to;
	edge_type et;

	directed_edge(pt::u32 f, pt::u32 t, edge_type et_)
	    : from(f), to(t), et(et_)
	{}
};

class BfsTree
{
	pt::u32 start; // index of the start vertex in vertices vector
	pt::u32 end;   // index of the end vertex in vertices vector

	std::vector<ptg::id_or_t> vertices;
	std::vector<directed_edge> edges_;

	// depth of each vertex in the BFS tree
	// idx is the BFS tree vertex index and value is the depth
	std::vector<pt::u32> depth;

	// ----
	// sort
	// ----
	// idx in the sort order are the sorts and values at that position are
	// indices in vertices vector
	std::vector<pt::u32> sort_order; // sort_order to BFS tree idx
	// idx in the tree to sort order
	std::vector<pt::u32> t_idx_to_sort_order;

public:
	// --------------
	// constructor(s)
	// --------------
	BfsTree() = default;

	// ---------
	// getter(s)
	// ---------

	[[nodiscard]] pt::u32 get_start() const
	{
		return this->start;
	}

	[[nodiscard]] pt::u32 get_end() const
	{
		return this->end;
	}

	[[nodiscard]] pt::u32 size() const
	{
		return static_cast<pt::u32>(this->vertices.size());
	}

	[[nodiscard]]
	std::vector<directed_edge> get_edges() const
	{
		return this->edges_;
	}

	[[nodiscard]]
	const ptg::id_or_t &get_vertex(pt::u32 t_v_idx) const
	{
		return this->vertices[t_v_idx];
	}

	[[nodiscard]] const std::vector<ptg::id_or_t> &get_vertices() const
	{
		return this->vertices;
	}

	[[nodiscard]] const std::vector<pt::u32> &get_sort_order() const
	{
		return this->sort_order;
	}

	[[nodiscard]]
	const ptg::id_or_t &get_sorted_vertex(pt::u32 sorted_idx) const
	{
		pt::u32 t_v_idx = this->sort_order[sorted_idx];
		return this->vertices[t_v_idx];
	}

	[[nodiscard]] pt::u32 get_sorted_pos(pt::u32 t_v_idx) const
	{
		if (this->t_idx_to_sort_order.size() <= t_v_idx)
			return pc::INVALID_IDX;

		return this->t_idx_to_sort_order[t_v_idx];
	}

	[[nodiscard]] pt::u32 get_depth(pt::u32 t_v_idx) const
	{
		return this->depth[t_v_idx];
	}

	[[nodiscard]] pt::u32 get_sort_idx(pt::u32 v_id) const
	{
		// i is t_v_idx
		for (pt::u32 i{}; i < this->sort_order.size(); i++) {
			auto [v_id_, _] = this->get_vertex(i);
			if (v_id_ == v_id)
				return this->get_sorted_pos(i);
		}

		return pc::INVALID_IDX;
	}

	[[nodiscard]]
	std::vector<pt::id_t> get_sorted() const
	{
		std::vector<pt::id_t> sorted;
		for (pt::u32 i{}; i < this->sort_order.size(); i++) {
			pt::u32 t_v_idx = this->sort_order[i];
			auto [v_id, _] = this->get_vertex(t_v_idx);
			sorted.push_back(v_id);
		}
		return sorted;
	}

	// ---------
	// setter(s)
	// ---------

	void set_start(pt::u32 s)
	{
		this->start = s;
	}

	void set_end(pt::u32 e)
	{
		this->end = e;
	}

	pt::u32 add_vertex(const ptg::id_or_t &v)
	{
		pt::u32 idx = this->vertices.size();
		this->vertices.push_back(v);
		return idx;
	}

	void add_tree_edge(pt::u32 from, pt::u32 to)
	{
		this->edges_.emplace_back(from, to, edge_type::tree_edge);
	}

	void add_cross_edge(pt::u32 from, pt::u32 to)
	{
		this->edges_.emplace_back(from, to, edge_type::cross_edge);
	}

	void print_dot(std::ostream &os) const
	{
		os << pv_cmp::format("digraph G {{\n"
				     "\trankdir = LR;\n"
				     "\tnode [shape = circle];\n"
				     "\tedge [arrowhead=vee];\n");

		// print vertices
		for (pt::u32 i{}; i < this->vertices.size(); i++) {
			os << pv_cmp::format(
				"\t{} [label=\"{} {} \n s:{} \n d:{}\" ", i, i,
				this->vertices[i].as_str(),
				this->get_sorted_pos(i), this->get_depth(i));

			if (i == this->start)
				os << "style=filled fillcolor=lightgreen";
			else if (i == this->end)
				os << "style=filled fillcolor=lightcoral";

			os << pv_cmp::format("];\n", i, i,
					     this->vertices[i].as_str(),
					     this->get_sorted_pos(i));
		}

		std::string edge_meta;
		for (const auto &[from, to, et] : this->edges_) {
			if (et == edge_type::tree_edge)
				edge_meta = R"( color="black" )";
			else if (et == edge_type::cross_edge)
				edge_meta = R"( style="dotted" color="red" )";

			os << pv_cmp::format("\t{} -> {} [{}];\n", from, to,
					     edge_meta);
		}

		os << "}\n";
	}

	// -----------
	// modifier(s)
	// -----------

	void comp_depth()
	{
		const pt::u32 N = this->size();
		this->depth = std::vector<pt::u32>(N, pc::INVALID_IDX);
		this->depth[this->start] = 0;

		std::queue<pt::u32> q;
		q.push(this->start);

		auto foo = [&](pt::u32 from, pt::u32 to,
			       pt::u32 curr_t_v_idx) -> bool
		{
			return from == curr_t_v_idx &&
			       this->depth[to] == pc::INVALID_IDX;
		};

		while (!q.empty()) {
			pt::u32 curr_t_v_idx = q.front();
			q.pop();
			pt::u32 curr_depth = this->depth[curr_t_v_idx];

			for (const auto &[from, to, et] : this->get_edges()) {
				if (et != edge_type::tree_edge)
					continue;

				if (foo(from, to, curr_t_v_idx)) {
					this->depth[to] = curr_depth + 1;
					q.push(to);
				}
			}
		}
	}

	// Modified Kahn's algorithm for
	// topological sort of the BFS tree from root to common
	// leaf using the given ref as a guide
	pt::status_t sort()
	{
		// bool dbg = s == ">7226>8008" ? true : false;

		const pt::u32 N = this->size();
		this->sort_order.reserve(N);
		this->t_idx_to_sort_order =
			std::vector<pt::u32>(N, pc::INVALID_IDX);

		std::vector<pt::u32> deg(N, 0);

		// compute in-degrees of all vertices
		for (const auto &[from, to, et] : this->get_edges()) {
			if (et == edge_type::tree_edge)
				deg[to]++;
			else if (et == edge_type::cross_edge && from != to)
				deg[to]++;
			// else if (et == edge_type::cross_edge && from != to &&
			//	 this->depth[from] < this->depth[to])
			//	deg[to]++;
		}

		std::queue<pt::u32> q;

		// start with vertices with in-degree 0, in this
		// case the root.
		q.push(this->start);

		std::set<pt::u32> done;
		done.insert(this->start);

		auto update_q = [&](pt::u32 u)
		{
			if (pv_cmp::contains(done, u))
				return;

			deg[u]--;

			if (deg[u] > 0)
				return;

			q.push(u);
			done.insert(u);
		};

		while (!q.empty()) {
			pt::u32 curr_t_v_idx = q.front();
			q.pop();
			pt::u32 sort_pos = this->sort_order.size();
			this->sort_order.emplace_back(curr_t_v_idx);
			this->t_idx_to_sort_order[curr_t_v_idx] = sort_pos;

			for (const auto &[from, to, et] : this->get_edges()) {
				if (et == edge_type::cross_edge && to == from)
					continue;

				if (from != curr_t_v_idx)
					continue;

				// if (et == edge_type::cross_edge &&
				//     this->depth[from] >= this->depth[to])
				//	continue;

				update_q(to);
			}
		}

		if (done.size() != N)
			return -1;

		return 0;
	}

	void set_sort(const std::list<pt::u32> &sorted_hap_idxs)
	{
		std::map<pt::u32, pt::u32> v_id_to_t_idx;
		for (pt::u32 i{}; i < this->size(); i++) {
			auto [v_id, _] = this->get_vertex(i);
			v_id_to_t_idx[v_id] = i;
		}

		const pt::u32 N = this->size();
		this->sort_order.resize(0);
		this->sort_order.reserve(N);
		this->t_idx_to_sort_order =
			std::vector<pt::u32>(N, pc::INVALID_IDX);

		pt::u32 i{};
		for (pt::u32 v_id : sorted_hap_idxs) {
			pt::u32 t_idx = v_id_to_t_idx[v_id];
			this->sort_order.emplace_back(t_idx);
			this->t_idx_to_sort_order[t_idx] = i;
			i++;
		}
	}
};

} // namespace ita::bfs

#endif // IT_BFS_TREE_HPP
