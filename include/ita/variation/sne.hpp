#ifndef IT_SNE_HPP
#define IT_SNE_HPP

#include <liteseq/refs.h> // for ref_walk, ref
#include <set>		  // for set
#include <utility>	  // for pair
#include <vector>	  // for vector

#include "ita/graph/slice_tree.hpp" // for it

#include "povu/common/constants.hpp"
#include "povu/common/core.hpp" // for pt, idx_t, id_t, op_t
#include "povu/graph/bidirected.hpp"
#include "povu/graph/types.hpp" // for or_e, id_or_t

namespace ita::sne // sne = seed and extend
{
constexpr std::string_view MODULE = "povu::overlay::sne";

namespace lq = liteseq;
namespace pgt = povu::types::graph;

struct pin {
	pt::u32 r_idx; // hap index
	pt::u32 idx;   // index in the hap walk
};

struct chain_link {
	pt::u32 idx;
	pt::u32 limit_left;
	pt::u32 limit_right;

	// hap_traversal dir; // direction of traversal

	chain_link static from_pin(const pin &p)
	{
		return chain_link{p.idx, pc::INVALID_IDX, pc::INVALID_IDX};
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
				pt::u32 N)
		{
			for (pt::u32 i = 0; i < N; i++) {
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
	// const pvr::RoV *rov_;
	std::vector<std::pair<pin, pin>> pin_pairs;
	std::map<pt::up_t<pt::u32>, std::set<pt::u32>> ref_to_pin_pairs;

public:
	// ---------------------
	// public constructor(s)
	// ---------------------
	// pin_cushion() = delete;

	explicit pin_cushion() = default;

	// ---------
	// getter(s)
	// ---------

	[[nodiscard]]
	bool is_empty() const
	{
		return this->pin_pairs.empty();
	}

	// [[nodiscard]]
	// std::string as_str() const
	// {
	//	return this->rov_->as_str();
	// }

	[[nodiscard]]
	std::vector<std::pair<pin, pin>>
	get_pin_pairs(pt::up_t<pt::u32> pp) const
	{
		if (!pv_cmp::contains(this->ref_to_pin_pairs, pp)) {
			return {};
		}

		std::vector<std::pair<pin, pin>> x;
		for (pt::u32 pos : this->ref_to_pin_pairs.at(pp)) {
			x.emplace_back(this->pin_pairs.at(pos));
		}

		return x;
	}

	// ---------
	// setter(s)
	// ---------

	void add_pin_pair(std::pair<pin, pin> &&pp)
	{
		pt::u32 pp_pos = this->pin_pairs.size();

		pt::u32 ref_h_idx = pp.first.r_idx;
		pt::u32 alt_h_idx = pp.second.r_idx;
		pt::up_t<pt::u32> pp_key{ref_h_idx, alt_h_idx};
		this->ref_to_pin_pairs[pp_key].insert(pp_pos);
		// this->ref_to_pin_pairs[alt_h_idx].insert(pp_pos);

		this->pin_pairs.emplace_back(pp);
	}
};

std::vector<ist::st> sne(const bd::VG &g, const pin_cushion &pcushions,
			 const std::set<pt::u32> &to_call_ref_ids);

}; // namespace ita::sne

// NOLINTNEXTLINE(misc-unused-alias-decls
namespace ise = ita::sne;

#endif // IT_SNE_HPP
