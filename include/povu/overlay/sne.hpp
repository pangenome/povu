#ifndef POVU_OVERLAY_SNE_HPP
#define POVU_OVERLAY_SNE_HPP

#include <liteseq/refs.h> // for ref_walk, ref
#include <optional>	  // for optional
#include <queue>
#include <set>	  // for set
#include <vector> // for vector

#include "povu/common/constants.hpp"
#include "povu/common/core.hpp" // for pt, idx_t, id_t, op_t
#include "povu/genomics/allele.hpp"
#include "povu/graph/bidirected.hpp"
#include "povu/graph/types.hpp"
#include "povu/variation/rov.hpp"

namespace povu::overlay::sne // sne = seed and extend
{
constexpr std::string_view MODULE = "povu::overlay::sne";

namespace lq = liteseq;
namespace pgt = povu::types::graph;

enum hap_traversal : pt::u8 {
	left,
	right
};

struct pin {
	pt::u32 r_idx;	   // ref index
	pt::u32 idx;	   // index in the ref
	hap_traversal dir; // direction of traversal
};

struct ref_sets {
	std::set<pt::u32> left_refs;
	std::set<pt::u32> right_refs;
};

struct pin_count {
	pt::u32 left_count = 0;
	pt::u32 right_count = 0;
};

struct chain_link {
	pt::u32 idx;
	pt::u32 limit_left;
	pt::u32 limit_right;
	hap_traversal dir; // direction of traversal

	chain_link static from_pin(const pin &p)
	{
		return chain_link{p.idx, pc::INVALID_IDX, pc::INVALID_IDX,
				  p.dir};
	}
};

// or matching overlap
struct extension {
	pt::u32 f_r_idx;
	pt::u32 r_r_idx;
	pt::u32 f_ref_start;
	pt::u32 r_ref_start;
	pt::u32 len;

	void dbg_print(const bd::VG &g) const
	{
		auto print = [](const lq::ref_walk *ref_w, pt::u32 ref_start,
				pt::u32 len)
		{
			for (pt::u32 i = 0; i < len; i++) {
				pt::idx_t ref_v_id =
					ref_w->v_ids[ref_start + i];
				pgt::or_e ref_o =
					ref_w->strands[ref_start + i] ==
							lq::strand::STRAND_FWD
						? pgt::or_e::forward
						: pgt::or_e::reverse;

				bd::id_or_t step{ref_v_id, ref_o};
				std::cerr << step.as_str();
			}
			std::cerr << "\n";
		};
		const lq::ref_walk *ref_w1 = g.get_ref_vec(f_r_idx)->walk;
		const lq::ref_walk *ref_w2 = g.get_ref_vec(r_r_idx)->walk;

		print(ref_w1, f_ref_start, len);
		print(ref_w2, r_ref_start - len + 1, len);
	}

	friend std::ostream &operator<<(std::ostream &os, const extension &ext)
	{
		os << "Extension(f_r_idx: " << ext.f_r_idx
		   << ", r_r_idx: " << ext.r_r_idx
		   << ", f_ref_start: " << ext.f_ref_start
		   << ", r_ref_start: " << ext.r_ref_start
		   << ", len: " << ext.len << ")";
		return os;
	}
};

struct link_pair {
	chain_link ref_link;
	chain_link alt_link;

	link_pair(chain_link rl, chain_link al) : ref_link(rl), alt_link(al)
	{}
};

struct chain_t {
	std::vector<link_pair> links;
	pt::u32 len_{};

	[[nodiscard]]
	pt::u32 len() const
	{
		return this->len_;
	}

	[[nodiscard]]
	bool empty() const
	{
		return this->len_ == 0;
	}
};

struct pin_cushion {
private:
	const pvr::RoV *rov_;
	// walk idx to pins in the walk
	std::map<pt::u32, std::vector<pin>> w_to_pins;

	// walk idx to refs that take the walk (left and right)
	std::map<pt::u32, ref_sets> w_to_refs;

	// walk to pin count left and pin count right
	std::map<pt::u32, pin_count> w_to_pin_counts;

	// ----------
	// updater(s)
	// ----------
	void update_w_to_pin_counts(pt::u32 w_idx, const pin &p)
	{
		auto &pc = this->w_to_pin_counts[w_idx];
		if (p.dir == hap_traversal::left)
			pc.left_count++;
		else
			pc.right_count++;
	}

	void update_w_to_refs(pt::u32 w_idx, const pin &p)
	{
		auto &[left_refs, right_refs] = this->w_to_refs[w_idx];
		if (p.dir == hap_traversal::left)
			left_refs.insert(p.r_idx);
		else
			right_refs.insert(p.r_idx);
	}

	void remove_walks(const std::set<pt::u32> &walks_to_remove)
	{
		for (pt::u32 w_idx : walks_to_remove) {
			this->w_to_pins.erase(w_idx);
			this->w_to_refs.erase(w_idx);
			this->w_to_pin_counts.erase(w_idx);
		}
	}

public:
	// ---------------------
	// public constructor(s)
	// ---------------------
	// pin_cushion() = default;

	explicit pin_cushion(const pvr::RoV *rov) : rov_(rov)
	{}

	// ---------
	// getter(s)
	// ---------
	[[nodiscard]]
	pt::u32 size() const
	{
		return static_cast<pt::u32>(this->w_to_pins.size());
	}

	[[nodiscard]]
	bool is_empty() const
	{
		return this->w_to_pins.empty();
	}

	[[nodiscard]]
	std::set<pt::u32> get_inv_walks() const
	{
		std::set<pt::u32> walks;
		for (const auto &[w_idx, pin_counts] : this->w_to_pin_counts) {
			auto [left_count, right_count] = pin_counts;
			if (left_count > 0 && right_count > 0)
				walks.insert(w_idx);
		}

		return walks;
	}

	[[nodiscard]]
	std::set<pt::u32> get_refs() const
	{
		std::set<pt::u32> refs;
		for (const auto &[w_idx, ref_set] : this->w_to_refs) {
			auto [left_refs, right_refs] = ref_set;
			refs.insert(left_refs.begin(), left_refs.end());
			refs.insert(right_refs.begin(), right_refs.end());
		}

		return refs;
	}

	[[nodiscard]]
	bool contains_ref(pt::u32 r_idx) const
	{
		for (const auto &[w_idx, ref_set] : this->w_to_refs) {
			auto [left_refs, right_refs] = ref_set;

			if (pv_cmp::contains(left_refs, r_idx) ||
			    pv_cmp::contains(right_refs, r_idx))
				return true;
		}

		return false;
	}

	[[nodiscard]]
	std::optional<hap_traversal> get_ref_direction(pt::u32 r_idx) const
	{
		for (const auto &[w_idx, ref_set] : this->w_to_refs) {
			auto [left_refs, right_refs] = ref_set;

			if (pv_cmp::contains(left_refs, r_idx))
				return hap_traversal::left;

			if (pv_cmp::contains(right_refs, r_idx))
				return hap_traversal::right;
		}

		return std::nullopt;
	}

	[[nodiscard]]
	std::optional<pin> get_ref_pin(pt::u32 r_idx) const
	{
		for (const auto &[_, pin_vec] : this->w_to_pins) {
			for (const pin &p : pin_vec) {
				if (p.r_idx == r_idx)
					return p;
			}
		}

		return std::nullopt;
	}

	[[nodiscard]]
	std::string as_str() const
	{
		return this->rov_->as_str();
	}

	// ---------
	// setter(s)
	// ---------

	void remove_uncallable_walks(const std::set<pt::u32> &to_call_ref_ids)
	{
		std::set<pt::u32> walks_to_remove;
		for (const auto &[w_idx, ref_set] : this->w_to_refs) {
			auto [left_refs, right_refs] = ref_set;

			bool has_callable_ref = false;
			for (const auto &r_idx : left_refs)
				if (pv_cmp::contains(to_call_ref_ids, r_idx))
					has_callable_ref = true;

			for (const auto &r_idx : right_refs)
				if (pv_cmp::contains(to_call_ref_ids, r_idx))
					has_callable_ref = true;

			// mark for removal if no callable refs
			if (!has_callable_ref)
				walks_to_remove.insert(w_idx);
		}

		this->remove_walks(walks_to_remove);
	}

	void add_pin(pt::u32 w_idx, pin &&p)
	{
		this->w_to_pins[w_idx].emplace_back(p);
		const pin &p_ = this->w_to_pins[w_idx].back();
		this->update_w_to_refs(w_idx, p_);
		this->update_w_to_pin_counts(w_idx, p_);
	}
};

struct edges {
	std::set<pt::u32> in;
	std::set<pt::u32> out;
};

class Tree
{
	std::map<int, std::vector<int>> adj; // Adjacency list representation
	std::map<int, int> child_to_parent_; // Parent map
	pt::u32 root;

public:
	// Add an edge (parent -> child)
	void add_edge(int parent, int child)
	{
		adj[parent].push_back(child);
		child_to_parent_[child] = parent;

		// Ensure child exists in adjacency list (even if it has no
		// children)
		if (adj.find(child) == adj.end())
			adj[child] = {};
	}

	[[nodiscard]]
	std::set<pt::u32> get_leaves() const
	{
		std::set<pt::u32> leaves;
		for (const auto &[node, children] : adj)
			if (children.empty())
				leaves.insert(node);

		return leaves;
	}

	[[nodiscard]]
	std::vector<pt::u32> get_stitch_path(pt::u32 leaf) const
	{
		// traverse until root or vertex with more than one child
		std::vector<pt::u32> path;
		pt::u32 current = leaf;

		while (true) {
			path.push_back(current);

			// check if current is root
			if (current == root)
				break;

			// get parent
			auto it = child_to_parent_.find(current);
			if (it == child_to_parent_.end())
				break; // no parent found, reached root

			pt::u32 parent = it->second;

			// check if parent has more than one child
			if (adj.at(parent).size() > 1)
				break;

			current = parent;
		}

		std::reverse(path.begin(), path.end());
		return path;
	}

	void set_root(pt::u32 r)
	{
		this->root = r;
	}

	// Perform BFS and print the tree level by level
	void bfs_traversal()
	{
		std::queue<pt::u32> q;		 // Queue for BFS
		std::map<pt::u32, bool> visited; // Track visited nodes

		q.push(this->root);
		visited[root] = true;

		while (!q.empty()) {
			pt::u32 node = q.front(); // Current node
			q.pop();

			std::cerr << node << " ";

			// Traverse all children of the current node
			for (pt::u32 child : adj[node]) {
				if (!visited[child]) {
					visited[child] = true;
					q.push(child);
				}
			}
		}

		std::cerr << "\n";
	}
};

class StitchGraph
{
	// key is vertex idx
	std::map<pt::u32, edges> e_;
	// std::map<pt::op_t<pt::u32>, pt::u32> v_to_v_idx_;
	std::vector<Tree> frst; // forest of BFS trees

	// Isolated nodes have no incoming or outgoing edges
	std::set<pt::u32> isolated;
	// Root nodes have no incoming edges
	std::set<pt::u32> roots;

public:
	void add_edge(pt::u32 from, pt::u32 to)
	{
		e_[from].out.insert(to);
		e_[to].in.insert(from);
	}

	void dbg_print_old()
	{
		// Print from each root using BFS
		for (pt::u32 root : this->roots) {
			std::cerr << "ROOT: " << root << "\n";
			std::queue<pt::u32> q;
			q.push(root);

			// Track visited nodes to avoid reprocessing
			// them
			std::set<pt::u32> visited;
			visited.insert(root);

			while (!q.empty()) {
				pt::u32 curr = q.front();
				q.pop();

				// Print all outgoing edges for the
				// current node
				const auto &[in, out] = e_.at(curr);
				for (pt::u32 o : out) {
					std::cerr << "    " << curr << " -> "
						  << o << "\n";
					if (!pv_cmp::contains(visited, o)) {
						q.push(o);
						visited.insert(o);
					}
				}
			}
		}
	}

	// vertices without incoming edges
	void comp_roots()
	{
		for (const auto &[v_idx, edges] : e_) {
			const auto &[in, out] = edges;
			if (in.empty() && !out.empty())
				this->roots.insert(v_idx);
			else if (in.empty() && out.empty())
				this->isolated.insert(v_idx);
		}
	}

	// Compute a BFS tree starting at a given root
	Tree comp_bfs_tree(pt::u32 start, std::set<pt::u32> &visited)
	{
		std::queue<pt::u32> q;
		Tree bfs_tree;

		q.push(start);
		visited.insert(start);
		bfs_tree.set_root(start);

		while (!q.empty()) {
			pt::u32 curr = q.front();
			q.pop();

			for (const auto &neighbor : e_[curr].out) {
				if (pv_cmp::contains(visited, neighbor))
					continue;

				visited.insert(neighbor);
				bfs_tree.add_edge(curr, neighbor);
				q.push(neighbor);
			}
		}

		return bfs_tree;
	}

	[[nodiscard]]
	const std::vector<Tree> &get_forest() const
	{
		return this->frst;
	}

	[[nodiscard]]
	const std::set<pt::u32> &get_isolated() const
	{
		return this->isolated;
	}

	// Compute BFS Forest (one tree per connected component)
	void compute_bfs_forest()
	{
		if (this->roots.empty())
			this->comp_roots();

		std::set<pt::u32> visited; // Track visited nodes
		for (pt::u32 r : this->roots)
			if (!pv_cmp::contains(visited, r))
				this->frst.push_back(comp_bfs_tree(r, visited));
	}

	void dbg_print()
	{
		pt::u32 tree_num = 1;
		for (auto &tree : this->frst) {
			std::cerr << "TREE " << tree_num++ << ":\n";
			tree.bfs_traversal(); // Assuming 0 is the root
		}
	}
};

void sne(const bd::VG &g, const std::vector<pin_cushion> &pcushions,
	 const std::set<pt::id_t> &to_call_ref_ids,
	 std::vector<pga::Exp> &exps);
}; // namespace povu::overlay::sne

// NOLINTNEXTLINE(misc-unused-alias-decls
namespace pos = povu::overlay::sne;

#endif // POVU_OVERLAY_SNE_HPP
