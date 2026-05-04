#ifndef PV_TREE_UTILS_HPP
#define PV_TREE_UTILS_HPP

#include <cstddef>     // for size_t
#include <iostream>    // for basic_ostream, operator<<, cerr
#include <map>	       // for map
#include <string>      // for operator<<
#include <string_view> // for string_view
#include <vector>      // for vector

#include <quilt/constants.hpp>	 // for pc
#include <quilt/graph_types.hpp> // for graph
#include <quilt/shim.hpp>	 // for format
#include <quilt/types.hpp>	 // for qt

#include "oza/graph/spanning_tree.hpp" // for Tree

namespace oza::tree_utils
{
inline constexpr std::string_view MODULE = "povu::tree_utils";

namespace pgt = quilt::types::graph;
namespace pst = oza::spanning_tree;

// key is edge idx of an edge and value is the v idx of the next braching vtx or
// of a leaf
using EdgeToBranch = std::map<qt::idx_t, qt::idx_t>;

struct BranchingMeta {
	// the black edge idx
	qt::idx_t black_e_idx = pc::INVALID_IDX;

	// branches sorted from the highest to lowest but if any child is black
	// then add first
	std::vector<qt::idx_t> sorted_br;

	// the vector should be sorted ...
	EdgeToBranch edge_data;
};

// key is a vertex idx and value is the branching meta data
using BranchDesc = std::map<qt::idx_t, BranchingMeta>;

struct tree_meta {
	std::vector<qt::idx_t> E;
	std::vector<qt::idx_t> D;

	std::vector<qt::idx_t>
		first; // idx is v_idx value is the first time it is seen in E
	std::vector<qt::idx_t> lo;  // LoA
	std::vector<qt::idx_t> HiD; // HiD

	std::map<qt::idx_t, qt::idx_t>
		pre; // idx is the pre-order the value is the v_idx
	std::map<qt::idx_t, qt::idx_t>
		post; // idx is the post-order the value is the v_idx

	// Gather all backedges
	std::vector<qt::idx_t> B;

	// prefix sum of the number of backedges
	std::vector<qt::idx_t> off;

	// a flat list of backedges
	std::vector<qt::idx_t> BE;

	std::vector<qt::idx_t> get_brackets(qt::idx_t v_idx) const
	{
		std::vector<qt::idx_t> brackets;
		qt::idx_t start = off[v_idx];
		qt::idx_t end = off[v_idx + 1];

		for (qt::idx_t i{start}; i < end; i++) {
			qt::idx_t be_idx = BE[i];
			brackets.push_back(be_idx);
		}

		return brackets;
	}

	std::vector<qt::idx_t> depth;

	void print()
	{

		// print lo
		std::cerr << "lo: \n";
		for (std::size_t v_idx = 0; v_idx < this->lo.size(); ++v_idx) {
			std::cerr << qs::format("({}, {}), ", v_idx,
						this->lo[v_idx]);
		}
		std::cerr << "\n\n";

		// print HiD
		std::cerr << "HiD: \n";
		for (std::size_t v_idx = 0; v_idx < this->HiD.size(); ++v_idx) {
			std::cerr << qs::format("({}, {}), ", v_idx,
						this->HiD[v_idx]);
		}

		// print E
		std::cerr << "E: \n";
		for (auto v_idx : this->E) {
			std::cerr << qs::format("{} ", v_idx);
		}
		std::cerr << "\n\n";

		// print D
		std::cerr << "D: \n";
		for (std::size_t v_idx = 0; v_idx < this->D.size(); ++v_idx) {
			std::cerr << qs::format("({}, {}), ", v_idx,
						this->D[v_idx]);
		}
		std::cerr << "\n\n";

		// print first
		std::cerr << "first (idx is v_idx value is the first time it "
			     "is seen in E): \n";
		for (qt::idx_t v_idx = 0; v_idx < this->first.size(); ++v_idx) {
			std::cerr << qs::format("({}, {}), ", v_idx,
						this->first[v_idx]);
		}
		std::cerr << "\n\n";

		// print depth
		std::cerr << "depth: \n";
		for (qt::idx_t v_idx = 0; v_idx < this->depth.size(); ++v_idx) {
			std::cerr << qs::format("({}, {}), ", v_idx,
						this->depth[v_idx]);
		}
	}
};

qt::idx_t find_lca(const tree_meta &tm, std::vector<qt::idx_t> &vtxs);
tree_meta gen_tree_meta(const pst::Tree &st);
BranchDesc br_desc(const pst::Tree &st);
} // namespace oza::tree_utils

#endif
