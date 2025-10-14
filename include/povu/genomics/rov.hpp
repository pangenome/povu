#ifndef POVU_GENOMICS_ROV_HPP
#define POVU_GENOMICS_ROV_HPP

#include <vector> // for vector

#include "povu/common/app.hpp" // for config
// #include "povu/common/constants.hpp"
#include "povu/common/core.hpp"
#include "povu/graph/bidirected.hpp" // for VG, bd
#include "povu/graph/pvst.hpp"	     // for Tree, VertexBase
#include "povu/graph/types.hpp"	     // for or_e, id_or_t, walk_t

namespace povu::genomics::rov
{
inline constexpr std::string_view MODULE = "povu::genomics::rov";
namespace pvst = povu::pvst;
namespace pgt = povu::types::graph;

enum class var_type_e : pt::u8 {
	del, // deletion
	ins, // insertion
	sub, // substitution
	und  // undetermined
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
	default:
		return "UND";
	}
}

std::ostream &operator<<(std::ostream &os, var_type_e vt);
var_type_e covariant(var_type_e a) noexcept;

struct raw_variant {
	pt::slice_t slice_a;
	pt::slice_t slice_b;
	var_type_e var_type;

	friend std::ostream &operator<<(std::ostream &os, const raw_variant &r)
	{
		os << " {" << r.slice_a << ", " << r.slice_b << ", "
		   << r.var_type << "}";
		return os;
	}

	friend bool operator==(const raw_variant &a, const raw_variant &b)
	{
		return a.slice_a == b.slice_a && a.slice_b == b.slice_b &&
		       a.var_type == b.var_type;
	}
};

struct pairwise_variants {
	pt::u32 walk_a;
	pt::u32 walk_b;
	std::vector<raw_variant> variants;

	[[nodiscard]]
	pt::u32 size() const noexcept
	{
		return static_cast<pt::u32>(this->variants.size());
	}

	[[nodiscard]]
	pairwise_variants(pt::u32 w1_idx, pt::u32 w2_idx) noexcept
	    : walk_a(w1_idx), walk_b(w2_idx), variants()
	{}

	void add_variant(raw_variant &&rv) noexcept
	{
		this->variants.emplace_back(rv);
	}

	friend std::ostream &operator<<(std::ostream &os,
					const pairwise_variants &pv)
	{
		os << "pairwise variants:\n";
		os << "(" << pv.walk_a << "," << pv.walk_b << ")\n";
		os << "variant count " << pv.variants.size() << "\n"
		   << "variants:\n";
		for (pt::u32 i{}; i < pv.variants.size(); i++) {
			os << "\t" << pv.variants[i];
			if (i != pv.variants.size() - 1)
				os << "\n";
		}

		return os;
	}
};

// struct walk_slice {
//	pt::u32 walk_idx;
//	pt::u32 start;
//	pt::u32 len;
// };

// // cannot be broken down further without losing information
// struct irreducible_rov {
//	std::vector<walk_slice> walk_slices;
//	// key is a pair of indexes in walk_slices, the value is the var type
//	std::map<pt::op_t<pt::u32>, var_type_e> var_types;
// };

/**
 * a collection of walks within a region of variation from start to end
 */
class RoV
{
	std::vector<pgt::walk_t> walks_;
	std::vector<pairwise_variants> irr_;
	const pvst::VertexBase *pvst_vtx;

public:
	// // print when RoV is moved
	// RoV(RoV &&other) noexcept
	//     : walks_(std::move(other.walks_)), irr_(std::move(other.irr_)),
	//       pvst_vtx(other.pvst_vtx)
	// {
	//	other.pvst_vtx = nullptr;
	// }

	// RoV(const RoV &other) = delete;		   // disable copy
	// constructor RoV &operator=(const RoV &other) = delete; // disable
	// copy assignment
	// // disable move assignment
	// RoV &operator=(RoV &&other) noexcept = delete;
	// ~RoV() = default;

	// --------------
	// constructor(s)
	// --------------

	RoV(const pvst::VertexBase *v) : walks_(), pvst_vtx(v)
	{}

	// ---------
	// getter(s)
	// ---------

	[[nodiscard]]
	pt::idx_t walk_count() const
	{
		return this->walks_.size();
	}

	[[nodiscard]]
	const pvst::VertexBase *get_pvst_vtx() const
	{
		return this->pvst_vtx;
	}

	[[nodiscard]]
	const std::vector<pgt::walk_t> &get_walks() const
	{
		return this->walks_;
	}

	[[nodiscard]]
	std::vector<pgt::walk_t> &get_walks_mut()
	{
		return this->walks_;
	}

	[[nodiscard]]
	const std::vector<pairwise_variants> &get_irreducibles() const
	{
		return this->irr_;
	}

	[[nodiscard]]
	const pgt::walk_t &get_walk(pt::idx_t idx) const
	{
		return this->walks_.at(idx);
	}

	// ---------
	// setter(s)
	// ---------

	void set_walks(std::vector<pgt::walk_t> &&walks)
	{
		this->walks_ = std::move(walks);
	}

	void add_irreducible(pairwise_variants &&irr_rov)
	{
		this->irr_.emplace_back(irr_rov);
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
 */
std::vector<RoV> gen_rov(const std::vector<pvst::Tree> &pvsts, const bd::VG &g,
			 const core::config &app_config);

#ifdef TESTING
void find_hidden(RoV &r);
#endif

} // namespace povu::genomics::rov

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pgr = povu::genomics::rov;

#endif
