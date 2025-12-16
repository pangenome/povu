#ifndef POVU_OVERLAY_INTERVAL_TREE_HPP
#define POVU_OVERLAY_INTERVAL_TREE_HPP

#include <liteseq/refs.h> // for ref_walk, ref
// #include <optional>	  // for optional
// #include <queue>
// #include <set> // for set
// #include <utility>
#include <map>
#include <vector> // for vector

#include "povu/common/constants.hpp"
#include "povu/common/core.hpp" // for pt, idx_t, id_t, op_t

// #include "povu/genomics/allele.hpp"
#include "povu/graph/bidirected.hpp"

// #include "povu/graph/types.hpp"
// #include "povu/variation/rov.hpp"

namespace povu::overlay::interval_tree
{
constexpr std::string_view MODULE = "povu::overlay::interval_tree";

namespace lq = liteseq;
namespace pgt = povu::types::graph;

struct alt {
	pt::u32 h_idx;
	pt::u32 h_start;
	pt::u32 len;
};

enum class update_type : pt::u8 {
	NO_OVERLAP = 0,
	EXISTS = 1,
	REPLACE = 2,
	MERGE = 3
};

enum class child_type : pt::u8 {
	LEFT = 0,
	RIGHT = 1
};

struct vertex {
	// hap data
	pt::u32 ref_h_start;

	// alt hap idx -> alt hap start
	std::map<pt::u32, std::vector<alt>> alts;

	// alt hap length -> haps that match the ref at this length
	std::map<pt::u32, std::set<pt::u32>> len_haps;

	// interval tree data
	// pt::u32 idx; // index of the vertex in the vertices array
	pt::u32 parent;
	pt::u32 left;
	pt::u32 right;

	static vertex create_root(pt::u32 ref_h_start, pt::u32 alt_h_idx,
				  pt::u32 alt_h_start, pt::u32 len)
	{
		vertex v;
		v.ref_h_start = ref_h_start;
		// v.idx = pc::INVALID_IDX;
		v.parent = pc::INVALID_IDX;
		v.left = pc::INVALID_IDX;
		v.right = pc::INVALID_IDX;

		alt a;
		a.h_idx = alt_h_idx;
		a.h_start = alt_h_start;
		a.len = len;
		v.alts[alt_h_idx].push_back(a);

		return v;
	}

	static vertex with_alt(pt::u32 p_idx, pt::u32 ref_h_start,
			       pt::u32 alt_h_idx, pt::u32 alt_h_start,
			       pt::u32 len)
	{
		vertex v;
		v.ref_h_start = ref_h_start;
		// v.idx = pc::INVALID_IDX;
		v.parent = p_idx;
		v.left = pc::INVALID_IDX;
		v.right = pc::INVALID_IDX;

		alt a;
		a.h_idx = alt_h_idx;
		a.h_start = alt_h_start;
		a.len = len;
		v.alts[alt_h_idx].push_back(a);

		return v;
	}

	[[nodiscard]]
	bool is_leaf() const
	{
		return this->left == pc::INVALID_IDX &&
		       this->right == pc::INVALID_IDX;
	}

	[[nodiscard]]
	std::map<pt::u32, std::vector<alt>> &get_alts_mut()
	{
		return this->alts;
	}

	[[nodiscard]]
	update_type has_overlap(pt::u32 alt_h_idx, pt::u32 alt_h_start,
				pt::u32 len)
	{
		if (!pv_cmp::contains(this->alts, alt_h_idx))
			return update_type::NO_OVERLAP;

		for (alt &a : this->alts[alt_h_idx]) {
			if (a.h_start == alt_h_start && a.len == len)
				return update_type::EXISTS; // already exists

			// what exists is contained
			if (alt_h_start < a.h_start &&
			    a.h_start + a.len <= alt_h_start + len) {
				// a.h_start = alt_h_start;
				// a.len = len;
				return update_type::REPLACE;
			}

			/* overlapping and need to merge*/

			if (a.h_idx + a.len >= alt_h_start &&
			    alt_h_start + len >= a.h_start) {
				return update_type::MERGE;
			}

			if (alt_h_start + len >= a.h_start &&
			    a.h_start + a.len >= alt_h_start) {
				return update_type::MERGE;
			}
		}
		return update_type::NO_OVERLAP;
	}

	void update_overlap(pt::u32 alt_h_idx, pt::u32 alt_h_start, pt::u32 len)
	{
		for (alt &a : this->alts[alt_h_idx]) {
			// already exists
			if (a.h_start == alt_h_start && a.len == len) {
				return;
			}

			// what is being added is contained
			if (a.h_start <= alt_h_start &&
			    alt_h_start + len <= a.h_start + a.len) {
				return;
			}

			// what exists is contained
			if (alt_h_start < a.h_start &&
			    a.h_start + a.len <= alt_h_start + len) {
				a.h_start = alt_h_start;
				a.len = len;
				return;
			}

			// overlapping, need to merge
			if (alt_h_start < a.h_start &&
			    a.h_idx + a.len >= alt_h_start) {
				// extend the existing alt
				a.h_start = alt_h_start;
				a.len = a.h_start - alt_h_start + a.len;
				return;
			}

			if (alt_h_start + len >= a.h_start) {
				// extend the existing alt
				a.h_start = alt_h_start; // new start
				a.len = a.h_start + a.len - alt_h_start;
				return;
			}
		}
	}

	void add_alt(pt::u32 alt_h_idx, pt::u32 alt_h_start, pt::u32 len)
	{
		alt a;
		a.h_idx = alt_h_idx;
		a.h_start = alt_h_start;
		a.len = len;

		this->alts[alt_h_idx].push_back(a);
		return;
	}

	[[nodiscard]]
	std::map<pt::u32, std::vector<alt>> get_same_len_alts() const
	{
		std::map<pt::u32, std::vector<alt>> alts_by_len;

		for (const auto &[alt_h_idx, alt_list] : this->alts) {
			for (const auto &a : alt_list) {
				alts_by_len[a.len].push_back(a);
			}
		}

		return alts_by_len;
	}

	void add_len_haps(pt::u32 h_idx, pt::u32 len)
	{
		this->len_haps[len].insert(h_idx);
	}

	[[nodiscard]]
	const std::set<pt::u32> &get_len_haps(pt::u32 len) const
	{
		if (!pv_cmp::contains(this->len_haps, len)) {
			ERR("No haps for length {}", len);
			std::exit(1);
		}

		return this->len_haps.at(len);
	}

	void set_left(pt::u32 v_idx)
	{
		if (this->left != pc::INVALID_IDX) {
			ERR("Left child already set for vertex {}",
			    this->ref_h_start);
			std::exit(1);
		}
		this->left = v_idx;
	}

	void clear_right()
	{
		this->right = pc::INVALID_IDX;
	}

	void clear_left()
	{
		this->left = pc::INVALID_IDX;
	}

	void set_right(pt::u32 v_idx)
	{
		if (this->right != pc::INVALID_IDX) {
			ERR("Right child already set for vertex {}",
			    this->ref_h_start);
			std::exit(1);
		}
		this->right = v_idx;
	}
};

struct it {
private:
	pt::u32 root;
	// std::vector<vertex> vertices;

	std::map<pt::u32, vertex> vertices;

	std::set<pt::u32> removed;

	pt::u32 ref_h_idx; // hap index
public:
	it() = delete;

	it(pt::u32 ref_h_idx_)
	    : root{pc::INVALID_IDX}, vertices{}, ref_h_idx{ref_h_idx_}
	{}

	[[nodiscard]]
	bool is_empty() const
	{
		return this->vertices.empty();
	}

	[[nodiscard]]
	pt::u32 get_ref_hap_idx() const
	{
		return this->ref_h_idx;
	}

	[[nodiscard]]
	pt::u32 get_root() const
	{
		return this->root;
	}

	[[nodiscard]]
	pt::u32 size() const
	{
		return this->vertices.size();
	}

	[[nodiscard]]
	const vertex &get_vertex(pt::u32 v_idx) const
	{
		return this->vertices.at(v_idx);
	}

	[[nodiscard]]
	const std::map<pt::u32, vertex> &get_vertices() const
	{
		return this->vertices;
	}

	[[nodiscard]]
	std::map<pt::u32, vertex> &get_vertices_mut()
	{
		return this->vertices;
	}

	[[nodiscard]]
	vertex &get_vertex_mut(pt::u32 v_idx)
	{
		return this->vertices.at(v_idx);
	}

	// modifiers

	void add_root(pt::u32 ref_h_start, pt::u32 alt_h_idx,
		      pt::u32 alt_h_start, pt::u32 len)
	{
		vertex v = vertex::create_root(ref_h_start, alt_h_idx,
					       alt_h_start, len);

		this->root = ref_h_start;
		this->vertices[ref_h_start] = v;
	}

	void add_non_root_leaf(pt::u32 p_idx, child_type ct,
			       pt::u32 ref_h_start, pt::u32 alt_h_idx,
			       pt::u32 alt_h_start, pt::u32 len)
	{
		vertex v = vertex::with_alt(p_idx, ref_h_start, alt_h_idx,
					    alt_h_start, len);

		this->vertices[ref_h_start] = v;

		/* update parent */
		switch (ct) {
		case child_type::LEFT:
			this->vertices[p_idx].set_left(ref_h_start);
			break;
		case child_type::RIGHT:
			this->vertices[p_idx].set_right(ref_h_start);
			break;
		}
	}

	void replace(pt::u32 curr_idx, pt::u32 ref_h_start, pt::u32 h2_idx,
		     pt::u32 h2_start, pt::u32 len)
	{
		vertex &curr_v = this->vertices[curr_idx];
		std::map<pt::u32, std::vector<alt>> &all_alts =
			curr_v.get_alts_mut();

		auto remove_overlap = [](std::vector<alt> &alts,
					 pt::u32 h2_start_, pt::u32 len_)
		{
			for (auto it_ = alts.begin(); it_ != alts.end();
			     ++it_) {
				// what is being added is contained

				if (h2_start_ < it_->h_start &&
				    it_->h_start + it_->len <=
					    h2_start_ + len_) {
					alts.erase(it_);
					break;
				}
			}
		};

		if (all_alts.size() > 1) {
			// std::cerr << "MULTI ALTS REPLACE\n";
			std::vector<alt> &alts = all_alts[h2_idx];
			remove_overlap(alts, h2_start, len);

			if (alts.empty())
				all_alts.erase(h2_idx);

			/* add new */
			this->add_vertex(ref_h_start, h2_idx, h2_start, len);
		}
		else { // size is 1
			// std::cerr << "SINGLE ALT REPLACE\n";
			std::vector<alt> &alts = all_alts[h2_idx];
			remove_overlap(alts, h2_start, len);

			if (alts.empty()) {
				this->remove_vertex(curr_idx);
				// std::cerr << "REMOVED VERTEX " << curr_idx
				//	  << " due to replace\n";
			}
			else {
				// std::cerr << "UPDATED EXISTING VERTEX "
				//	  << curr_idx
				//	  << " due to "
				//	     "replace\n";
			}

			this->add_vertex(ref_h_start, h2_idx, h2_start, len);
		}
	}

	void remove_root()
	{
		if (this->root == pc::INVALID_IDX)
			return;

		vertex &root_v = this->get_vertex_mut(this->root);
		if (root_v.left != pc::INVALID_IDX) {
			this->root = root_v.left;
			vertex &new_root = this->get_vertex_mut(this->root);
			new_root.parent = pc::INVALID_IDX;

			if (root_v.right != pc::INVALID_IDX) {
				// attach right subtree
				pt::u32 current_idx = this->root;
				vertex &current_v =
					this->get_vertex_mut(current_idx);
				while (current_v.right != pc::INVALID_IDX) {
					current_idx = current_v.right;
					current_v = this->get_vertex_mut(
						current_idx);
				}
				current_v.set_right(root_v.right);
				vertex &right_v =
					this->get_vertex_mut(root_v.right);
				right_v.parent = current_idx;
			}
		}
		else if (root_v.right != pc::INVALID_IDX) {
			this->root = root_v.right;
			vertex &new_root = this->get_vertex_mut(this->root);
			new_root.parent = pc::INVALID_IDX;
		}
		else {
			// tree is now empty
			this->root = pc::INVALID_IDX;
		}

		this->vertices.erase(root_v.ref_h_start);
	}

	void remove_leaf(pt::u32 v_idx)
	{
		// std::cerr << "REMOVING LEAF VERTEX " << v_idx << "\n";

		if (this->root == v_idx) { // the leaf is the root
			this->remove_root();
			return;
		}

		// update parent to remove reference to this leaf
		pt::u32 p_idx = this->get_vertex(v_idx).parent;
		vertex &p_v = this->get_vertex_mut(p_idx);
		if (p_v.left == v_idx)
			p_v.clear_left();
		else if (p_v.right == v_idx)
			p_v.clear_right();
		else {
			ERR("Inconsistent tree state when removing leaf");
			std::exit(1);
		}

		this->vertices.erase(v_idx);
	}

	void remove_internal(pt::u32 v_idx)
	{

		// std::cerr << "REMOVING INTERNAL VERTEX " << v_idx << "\n";

		vertex &v = this->get_vertex_mut(v_idx);

		if (v.left != pc::INVALID_IDX && v.right != pc::INVALID_IDX) {
			// two children
			// find inorder successor (leftmost of right subtree)
			pt::u32 succ_idx = v.right;
			vertex &succ_v = this->get_vertex_mut(succ_idx);
			while (succ_v.left != pc::INVALID_IDX) {
				succ_idx = succ_v.left;
				succ_v = this->get_vertex_mut(succ_idx);
			}

			// copy successor data to current vertex
			v.ref_h_start = succ_v.ref_h_start;
			v.alts = succ_v.alts;
			v.len_haps = succ_v.len_haps;

			// remove successor
			if (succ_v.is_leaf()) {
				remove_leaf(succ_idx);
			}
			else {
				remove_internal(succ_idx);
			}
		}
		else if (v.left != pc::INVALID_IDX ||
			 v.right != pc::INVALID_IDX) {
			// one child
			pt::u32 child_idx =
				(v.left != pc::INVALID_IDX) ? v.left : v.right;
			vertex &child_v = this->get_vertex_mut(child_idx);

			// update parent to point to child
			pt::u32 p_idx = v.parent;
			if (p_idx != pc::INVALID_IDX) {
				vertex &p_v = this->get_vertex_mut(p_idx);
				if (p_v.left == v_idx)
					p_v.left = child_idx;
				else if (p_v.right == v_idx)
					p_v.right = child_idx;
			}
			child_v.parent = p_idx;

			// this->removed.insert(v_idx);
		}
		else {
			ERR("Inconsistent tree state when removing internal");
			std::exit(1);
		}
	}

	void remove_vertex(pt::u32 v_idx)
	{
		if (this->get_vertex(v_idx).is_leaf())
			remove_leaf(v_idx);
		else
			remove_internal(v_idx);
	}

	void add_non_root(pt::u32 ref_h_start, pt::u32 alt_h_idx,
			  pt::u32 alt_h_start, pt::u32 len)
	{
		// find a matching ref_h_start in the tree
		pt::u32 current_idx = this->root;
		while (current_idx != pc::INVALID_IDX) {
			vertex &current = this->vertices[current_idx];

			update_type ut = current.has_overlap(alt_h_idx,
							     alt_h_start, len);

			// if (ut != update_type::NO_OVERLAP) {
			//	std::cerr << "UPDATING OVERLAP at curr "
			//		     "ref_h_start "
			//		  << current.ref_h_start
			//		  << " new ref h start " << ref_h_start
			//		  << " for alt_h_idx " << alt_h_idx
			//		  << ", alt_h_start " << alt_h_start
			//		  << ", len " << len << "\n";
			// }

			switch (ut) {
			case update_type::REPLACE:
				// std::cerr << "Replace\n";
				this->replace(current_idx, ref_h_start,
					      alt_h_idx, alt_h_start, len);
				return;
			case update_type::MERGE:
				// std::cerr << "Merge\n";
				current.update_overlap(alt_h_idx, alt_h_start,
						       len);
				return;
			case update_type::EXISTS:
				// std::cerr << "Exists\n";
				return;
			case update_type::NO_OVERLAP: // continue
				break;
			}

			if (ref_h_start < current.ref_h_start) {
				// go left
				if (current.left == pc::INVALID_IDX) {
					// std::cerr << "INSERTING at LEFT of "
					//	     "ref_h_start "
					//	  << current.ref_h_start
					//	  << "\n";
					// insert here
					this->add_non_root_leaf(
						current_idx, child_type::LEFT,
						ref_h_start, alt_h_idx,
						alt_h_start, len);
					return;
				}
				else {
					current_idx = current.left;
				}
			}
			else if (ref_h_start > current.ref_h_start) {
				// go right
				if (current.right == pc::INVALID_IDX) {
					// std::cerr << "INSERTING at RIGHT of "
					//	     "ref_h_start "
					//	  << current.ref_h_start
					//	  << "\n";
					// insert here
					this->add_non_root_leaf(
						current_idx, child_type::RIGHT,
						ref_h_start, alt_h_idx,
						alt_h_start, len);
					return;
				}
				else {
					current_idx = current.right;
				}
			}
			else {
				// std::cerr << "EQ\n";
				// matching ref_h_start found, add alt
				current.add_alt(alt_h_idx, alt_h_start, len);
				return;
			}
		}
	}

	void add_vertex(pt::u32 ref_h_start, pt::u32 alt_h_idx,
			pt::u32 alt_h_start, pt::u32 len)

	{
		if (this->root == pc::INVALID_IDX) {
			this->root = this->vertices.size();
			this->add_root(ref_h_start, alt_h_idx, alt_h_start,
				       len);
			return;
		}

		this->add_non_root(ref_h_start, alt_h_idx, alt_h_start, len);
	}

	void dbg_print_vertex(std::ostream &os, const bd::VG &g,
			      const vertex &v) const
	{
		auto print_slice =
			[&](const lq::ref_walk *h_w, pt::u32 start, pt::u32 len)
		{
			pt::u32 N = start + len;
			for (pt::u32 i = start; i < N; i++) {
				pt::idx_t ref_v_id = h_w->v_ids[i];
				pgt::or_e ref_o =
					h_w->strands[i] ==
							lq::strand::STRAND_FWD
						? pgt::or_e::forward
						: pgt::or_e::reverse;
				bd::id_or_t step{ref_v_id, ref_o};
				std::cerr << step.as_str();
			}
			os << "\n";
		};

		for (const auto &[len, alts] : v.get_same_len_alts()) {
			const lq::ref_walk *ref_h_w =
				g.get_ref_vec(ref_h_idx)->walk;
			os << "REF: ";
			print_slice(ref_h_w, v.ref_h_start, len);
			for (const auto &a : alts) {
				const lq::ref_walk *alt_h_w =
					g.get_ref_vec(a.h_idx)->walk;
				os << "  ALT (h_idx: " << a.h_idx << "): ";
				print_slice(alt_h_w, a.h_start, len);
			}
		}
	}

	void dbg_print(std::ostream &os, const bd::VG &g) const
	{
		auto print_slice =
			[&](const lq::ref_walk *h_w, pt::u32 start, pt::u32 len)
		{
			pt::u32 N = start + len;
			for (pt::u32 i = start; i < N; i++) {
				pt::idx_t ref_v_id = h_w->v_ids[i];
				pgt::or_e ref_o =
					h_w->strands[i] ==
							lq::strand::STRAND_FWD
						? pgt::or_e::forward
						: pgt::or_e::reverse;
				bd::id_or_t step{ref_v_id, ref_o};
				std::cerr << step.as_str();
			}
			os << "\n";
		};

		os << "\n-------------------------\n";
		os << "Interval Tree (ref_h_idx: " << this->ref_h_idx << ")\n";
		for (const auto &[_, v] : this->vertices) {

			for (const auto &[len, alts] : v.get_same_len_alts()) {
				const lq::ref_walk *ref_h_w =
					g.get_ref_vec(ref_h_idx)->walk;
				os << "REF: ";
				print_slice(ref_h_w, v.ref_h_start, len);
				for (const auto &a : alts) {
					const lq::ref_walk *alt_h_w =
						g.get_ref_vec(a.h_idx)->walk;
					os << "  ALT (h_idx: " << a.h_idx
					   << "): ";
					print_slice(alt_h_w, a.h_start, len);
				}
			}
		}
		os << "\n-------------------------\n";
	}
};

} // namespace povu::overlay::interval_tree

namespace poi = povu::overlay::interval_tree;

#endif // POVU_OVERLAY_INTERVAL_TREE_HPP
