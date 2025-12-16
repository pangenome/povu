#ifndef POVU_ROV_HPP
#define POVU_ROV_HPP

#include <optional> // for optional
#include <string>   // for string
#include <vector>   // for vector

#include "povu/common/core.hpp"	     // for pt
#include "povu/graph/bidirected.hpp" // for VG, bd
#include "povu/graph/pvst.hpp"	     // for Tree, VertexBase
#include "povu/graph/types.hpp"	     // for or_e, id_or_t, walk_t

namespace povu::var::rov
{
inline constexpr std::string_view MODULE = "povu::genomics::rov";
namespace pvst = povu::pvst;
namespace pgt = povu::types::graph;

enum class var_type_e : pt::u8 {
	del, // deletion
	ins, // insertion
	sub, // substitution
	inv, // inversion
};

constexpr std::string_view to_string_view(var_type_e vt) noexcept
{
	switch (vt) {
	case var_type_e::del:
		return "DEL";
	case var_type_e::ins:
		return "INS";
	case var_type_e::sub:
		return "SUB";
	case var_type_e::inv:
		return "INV";
	}

	ERR("Unknown variant type");

	return "UNKNOWN";
}

std::ostream &operator<<(std::ostream &os, var_type_e vt);
var_type_e covariant(var_type_e a) noexcept;

/**
 * Represents a genomic region for filtering RoVs
 * Format: ref_name:start-end (e.g., "chr17:43044294-43125482")
 */
struct genomic_region {
	std::string ref_name; // reference path name
	pt::idx_t start;      // start position (0-based)
	pt::idx_t end;	      // end position (exclusive)

	genomic_region() : ref_name(""), start(0), end(0)
	{}

	genomic_region(const std::string &name, pt::idx_t s, pt::idx_t e)
	    : ref_name(name), start(s), end(e)
	{}

	[[nodiscard]]
	bool is_valid() const
	{
		return !ref_name.empty() && start < end;
	}
};

/**
 * Parse a genomic region string (format: ref:start-end)
 * @param region_str The region string to parse
 * @return Optional genomic_region if parsing succeeds, std::nullopt otherwise
 */
std::optional<genomic_region>
parse_genomic_region(const std::string &region_str);

struct enhanced_walk {
	pgt::walk_t steps;
	std::vector<pt::op_t<pt::u32>> cycles;
};

/**
 * a collection of walks within a region of variation from start to end
 */
class RoV
{
	// TODO: make this a variant
	// -------------------------
	//
	// (a) a set of walks enumerated
	// std::vector<enhanced_walk> e_walks_;
	// (b) a BFS tree
	// povu::graph::bfs::BfsTree bfs_tree;

	std::vector<pt::id_t> vertices;
	// v_id to idx in vertices (sort order)
	std::map<pt::u32, pt::u32> sort_order;

	const pvst::VertexBase *pvst_vtx;

public:
	// --------------
	// constructor(s)
	// --------------

	RoV(const pvst::VertexBase *v) : pvst_vtx(v)
	{}

	// ---------
	// getter(s)
	// ---------

	[[nodiscard]]
	const pvst::VertexBase *get_pvst_vtx() const
	{
		return this->pvst_vtx;
	}

	// [[nodiscard]]
	// const povu::graph::bfs::BfsTree &get_bfs_tree() const
	// {
	//	return this->bfs_tree;
	// }

	// [[nodiscard]]
	// povu::graph::bfs::BfsTree &get_bfs_tree_mut()
	// {
	//	return this->bfs_tree;
	// }

	[[nodiscard]]
	pt::u32 size() const noexcept
	{
		return static_cast<pt::u32>(this->vertices.size());
	}

	[[nodiscard]]
	pt::u32 get_vertex_count() const noexcept
	{
		return static_cast<pt::u32>(this->vertices.size());
	}

	[[nodiscard]]
	bool can_be_non_planar() const
	{
		pvst::vf_e fam = this->pvst_vtx->get_fam();
		if (fam == pvst::vf_e::tiny || fam == pvst::vf_e::parallel)
			return false;

		return true;
	}

	[[nodiscard]]
	bool is_sorted() const noexcept
	{
		return !this->vertices.empty();
	}

	[[nodiscard]]
	pt::u32 get_sorted_vertex(pt::u32 i) const
	{
		return this->vertices.at(i);
	}

	[[nodiscard]]
	const std::vector<pt::id_t> &get_sorted_vertices() const
	{
		return this->vertices;
	}

	[[nodiscard]]
	pt::u32 get_sorted_pos(pt::id_t v_id) const
	{
		// std::cerr << "Getting sorted pos for v_id: " << v_id << "\n";
		if (!pv_cmp::contains(this->sort_order, v_id))
			return pc::INVALID_IDX;

		return this->sort_order.at(v_id);
	}

	// ---------
	// setter(s)
	// ---------
	// takes an iterable
	template <typename InputIt>
	void add_sort_data(InputIt first, InputIt last)
	{
		for (InputIt it = first; it != last; ++it) {
			pt::id_t v_id = *it;
			this->sort_order[v_id] = this->vertices.size();
			this->vertices.push_back(v_id);
		}
	}

	// --------
	// other(s)
	// --------

	[[nodiscard]]
	std::string as_str() const
	{
		return this->pvst_vtx->as_str();
	}
};

/**
 * find walks in the graph based on the leaves of the pvst
 * initialize RoVs from flubbles
 * @param region Optional genomic region to filter RoVs
 */
std::vector<RoV>
gen_rov(const std::vector<pvst::Tree> &pvsts, const bd::VG &g,
	const std::set<pt::id_t> &to_call_ref_ids,
	const std::optional<genomic_region> &region = std::nullopt);

} // namespace povu::var::rov

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pvr = povu::var::rov;

#endif // POVU_ROV_HPP
