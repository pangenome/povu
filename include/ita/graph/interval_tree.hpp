#ifndef IT_IT_TREE_HPP
#define IT_IT_TREE_HPP

#include <cstdlib>
#include <map>
#include <optional>
#include <set>

#include "povu/common/constants.hpp"
#include "povu/common/core.hpp" // for pt, idx_t, id_t, op_t
#include "povu/common/log.hpp"	// for ERR

namespace ita::interval_tree
{
constexpr std::string_view MODULE = "ita::interval_tree";

/* insert specific */
enum class child_type : pt::u8 {
	LEFT = 0,
	RIGHT = 1
};

struct vertex {
	pt::u32 fwd_start;
	pt::u32 rev_start;
	pt::u32 len;

	pt::u32 left_child{pc::INVALID_IDX};
	pt::u32 right_child{pc::INVALID_IDX};

	vertex() = delete;

	vertex(pt::u32 fwd_s, pt::u32 rev_s, pt::u32 l)
	    : fwd_start{fwd_s}, rev_start{rev_s}, len{l}
	{}

	void add_child(pt::u32 child_idx, child_type ct)
	{
		if (ct == child_type::LEFT)
			this->left_child = child_idx;
		else
			this->right_child = child_idx;
	}

	[[nodiscard]]
	bool has_child_left() const
	{
		return this->left_child != pc::INVALID_IDX;
	}

	[[nodiscard]]
	bool has_child_right() const
	{
		return this->right_child != pc::INVALID_IDX;
	}
};

struct interval_tree {
private:
	pt::u32 ref_h_idx_{pc::INVALID_IDX}; // ref hap index
	pt::u32 root_{pc::INVALID_IDX};
	std::map<pt::u32, vertex> vertices{};

public:
	interval_tree() = delete;

	interval_tree(pt::u32 hap_idx) : ref_h_idx_{hap_idx}
	{}

	[[nodiscard]]
	bool empty() const
	{
		return this->vertices.empty();
	}

	[[nodiscard]]
	pt::u32 get_ref_hap_idx() const
	{
		return this->ref_h_idx_;
	}

	[[nodiscard]]
	const std::map<pt::u32, vertex> &get_intervals() const
	{
		return this->vertices;
	}

	void add_root(pt::u32 fwd_start, pt::u32 rev_start, pt::u32 len)
	{
		vertex v{fwd_start, rev_start, len};
		this->vertices.emplace(fwd_start, v);
		this->root_ = fwd_start;
	}

	void add_internal(pt::u32 fwd_start, pt::u32 rev_start, pt::u32 len)
	{
		pt::u32 curr_idx = this->root_;
		while (true) {
			vertex &curr_v = this->vertices.at(curr_idx);

			// extend current vertex
			// overlap
			if (curr_v.fwd_start + curr_v.len >= fwd_start &&
			    fwd_start + len > curr_v.fwd_start + curr_v.len) {
				pt::u32 new_end = fwd_start + len;
				pt::u32 new_len = new_end - curr_v.fwd_start;
				curr_v.len = std::max(curr_v.len, new_len);
				return;
			}

			if (fwd_start < curr_v.fwd_start) {
				if (curr_v.has_child_left()) {
					curr_idx = curr_v.left_child;
				}
				else {
					vertex new_v{fwd_start, rev_start, len};
					this->vertices.emplace(fwd_start,
							       new_v);
					curr_v.add_child(fwd_start,
							 child_type::LEFT);
					return;
				}
			}
			else { // fwd_start >= curr_v.fwd_start
				if (curr_v.has_child_right()) {
					curr_idx = curr_v.right_child;
				}
				else {
					vertex new_v{fwd_start, rev_start, len};
					this->vertices.emplace(fwd_start,
							       new_v);
					curr_v.add_child(fwd_start,
							 child_type::RIGHT);
					return;
				}
			}
		}
	}

	void add(pt::u32 fwd_start, pt::u32 rev_start, pt::u32 len)
	{
		if (this->empty())
			this->add_root(fwd_start, rev_start, len);
		else
			this->add_internal(fwd_start, rev_start, len);
	}
};

} // namespace ita::interval_tree

namespace iit = ita::interval_tree;

#endif // IT_IT_TREE_HPP
