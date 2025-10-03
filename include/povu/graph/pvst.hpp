#ifndef POVU_TYPES_PVST_HPP
#define POVU_TYPES_PVST_HPP

#include <algorithm>
#include <optional>
#include <stack>
#include <string_view>

"#include "povu/common/compat.hpp"
"#include "povu/common/constants.hpp"
"#include "povu/common/log.hpp"
#include "types.hpp"

/* === PVST pangenome variation structure tree === */

namespace povu::pvst
{
inline constexpr std::string_view MODULE = "povu::pvst";

namespace pc = povu::constants;
namespace pgt = povu::types::graph;

/*
hierarchy of the types of vertices
vertex clan
  ↳ vertex family

the dummy doesn't belong to any clan or family

vertex family:
  flubble
    ↳ tiny
    ↳ parallel
    ↳ generic flubble
  concealed
    ↳ smothered
  midi bubble

vertex clan:
fl_like
  ↳ flubble
  ↳ tiny
  ↳ parallel
subflubble
  ↳ concealed
  ↳ smothered
  ↳ midi bubble
*/

// prototype the clan enum
enum class vertex_clan_e;
typedef vertex_clan_e vc_e;

enum class vertex_family_e {
	dummy,
	flubble, // TODO: aka generic flubble, rename to generic_flubble?
	tiny,
	parallel,
	concealed,
	smothered,
	midi,
};
typedef vertex_family_e vt_e;
typedef vertex_family_e vf_e;

enum class vertex_clan_e {
	fl_like, // flubble-like vertices
	subflubble,
};

inline const char *to_str(vf_e f)
{
	switch (f) {
	case vf_e::dummy:
		return "dummy";
	case vf_e::flubble:
		return "flubble";
	case vf_e::tiny:
		return "tiny";
	case vf_e::parallel:
		return "parallel";
	case vf_e::concealed:
		return "concealed";
	case vf_e::smothered:
		return "smothered";
	case vf_e::midi:
		return "midi";
	default:
		return "unknown";
	}
}

inline std::ostream &operator<<(std::ostream &os, vf_e t)
{
	return os << to_str(t);
}

/**
 * Deterimines the vertex clan based on its family.
 * @param t vertex family
 * @return an optional vertex clan
 *
 * This function is used to determine the category of a vertex based on its
 * type. It is useful for categorizing vertices in the PVST structure.
 *
 * @note The function returns a nullopt if the vertex family does not belong to
 * any defined clan.
 */
[[nodiscard]] inline std::optional<vc_e> to_clan(vf_e f) noexcept
{
	switch (f) {
	case vf_e::flubble:
	case vf_e::tiny:
	case vf_e::parallel:
		return vc_e::fl_like;
	case vf_e::concealed:
	case vf_e::smothered:
	case vf_e::midi:
		return vc_e::subflubble;
	default:
		return std::nullopt; // undefined category
	}
}

struct endpoints {
	pgt::id_or_t s;
	pgt::id_or_t e;
};

// ==========================
// Related to graph traversal
// ==========================

// direction of the traversal
enum class route_e {
	s2e, // start to end
	e2s, // end to start
};
typedef route_e rt_e;

inline const char *to_str(route_e r)
{
	switch (r) {
	case route_e::s2e:
		return "L";
	case route_e::e2s:
		return "R";
	default:
		ERR("Unknown route_e value: {}", static_cast<int>(r));
		std::exit(1);
	}
}

inline std::ostream &operator<<(std::ostream &os, route_e r)
{
	return os << to_str(r);
}

// bounds in the bidirected graph
struct route_params_t {
	pgt::id_or_t start;
	pgt::id_or_t end;
	rt_e route;
};

typedef route_params_t rp_t;

// ============================
// fouces on the flubble within the ST
// ==========================

// where the concealed vertex is located in relation to its parent flubble
enum class concealed_fl_location_e {
	ai_trunk,
	ai_branch,
	zi_trunk,
	zi_branch,
	undefined
};
typedef concealed_fl_location_e cl_e;

// type of the concealed vertex
enum class concealed_fl_boundary_e {
	g,
	s
};
typedef concealed_fl_boundary_e cb_e;

// bounds in the spanning tree
struct bounds_t {
	pt::idx_t upper;
	pt::idx_t lower; // when invalid all leaves are the lower boundaries
};

const bounds_t INVALID_BOUNDS{pc::INVALID_IDX, pc::INVALID_IDX};

/* an abstract class for vertices  */
class VertexBase
{
	pt::id_t idx_;	   // idx of the vertex in the pvst
	pt::idx_t height_; // height of the vertex in the pvst
	vf_e fam_;

public:
	// --------------
	// constructor(s)
	// --------------
	VertexBase(pt::id_t idx, vf_e fam)
	    : idx_(idx), height_(pc::INVALID_IDX), fam_(fam)
	{}

	// -------
	// getters
	// -------
	pt::id_t get_idx() const
	{
		return this->idx_;
	}

	vf_e get_fam() const
	{
		return this->fam_;
	}

	pt::idx_t get_height() const
	{
		return this->height_;
	}

	// ----------------------
	// pure virtual functions
	// ----------------------
	virtual std::string as_str() const = 0;
	// TODO: make this result a const reference?
	virtual std::optional<route_params_t> get_route_params() const = 0;
	virtual ~VertexBase() = default;

	// -------
	// setters
	// -------
	void set_idx(pt::id_t idx)
	{
		this->idx_ = idx;
	}

	void set_type(vf_e type)
	{
		this->fam_ = type;
	}

	void set_height(pt::idx_t height)
	{
		this->height_ = height;
	}
};

class Dummy : public VertexBase
{
public:
	// --------------
	// constructor(s)
	// --------------
	Dummy() : VertexBase(pc::INVALID_IDX, vf_e::dummy)
	{}

	// -------
	// getters
	// -------
	std::string as_str() const override
	{
		return ".";
	}

	std::optional<route_params_t> get_route_params() const override
	{
		return std::nullopt; // dummy vertex does not have route params
	}
};

class Flubble : public VertexBase
{
	pgt::id_or_t a_; // start
	pgt::id_or_t z_; // end

	pt::idx_t ai_;
	pt::idx_t zi_;

	pt::idx_t m_;
	pt::idx_t n_;

	// a flubble route is always L to R
	const route_e route_{route_e::s2e};

	// --------------
	// constructor(s)
	// --------------
	Flubble(vf_e vf, pgt::id_or_t a, pgt::id_or_t z, pt::idx_t ai,
		pt::idx_t zi)
	    : VertexBase(pc::INVALID_IDX, vf), a_(a), z_(z), ai_(ai), zi_(zi),
	      m_(pc::INVALID_IDX), n_(pc::INVALID_IDX)
	{}

public:
	// --------------
	// factory fns
	// --------------

	static Flubble create(pgt::id_or_t a, pgt::id_or_t z, pt::idx_t ai,
			      pt::idx_t zi)
	{
		return Flubble(vf_e::flubble, a, z, ai, zi);
	}

	static Flubble parse(vf_e f, const route_params_t &rp)
	{
		auto [a, z, _] = rp;
		pt::idx_t ai{pc::INVALID_IDX};
		pt::idx_t zi{pc::INVALID_IDX};
		return Flubble(f, a, z, ai, zi);
	}

	// -------
	// getters
	// -------
	pgt::id_or_t get_a() const
	{
		return this->a_;
	}

	pgt::id_or_t get_z() const
	{
		return this->z_;
	}

	pt::idx_t get_ai() const
	{
		return this->ai_;
	}

	pt::idx_t get_zi() const
	{
		return this->zi_;
	}

	pt::idx_t get_m() const
	{
		return this->m_;
	}

	pt::idx_t get_n() const
	{
		return this->n_;
	}

	bounds_t get_bounds() const
	{
		return {this->get_ai(), this->get_zi()};
	}

	std::optional<route_params_t> get_route_params() const override
	{
		return route_params_t{this->get_a(), this->get_z(),
				      this->route_};
	}

	// -------
	// setters
	// -------
	void set_m(pt::idx_t m)
	{
		this->m_ = m;
	}

	void set_n(pt::idx_t n)
	{
		this->n_ = n;
	}

	// ------
	// others
	// ------
	std::string as_str() const override
	{
		return pv_cmp::format("{}{}", this->a_.as_str(),
				      this->z_.as_str());
	}
};

class Concealed : public VertexBase
{
	// TODO: rename to parent idx
	pt::idx_t fl_idx_;  // v idx of the parent flubble in the PVST
	cl_e location_;	    // type of the slubble (trunk or branch)
	pt::idx_t loc_idx_; // idx in the spanning tree for slubble

	// b for boundary
	pgt::id_or_t fl_b_; // a or z
	pgt::id_or_t cn_b_; // g or s

	bounds_t bounds_;

	route_e route_;

private:
	bool with_ai() const
	{
		return this->location_ == cl_e::ai_trunk ||
		       this->location_ == cl_e::ai_branch;
	}

	bool with_zi() const
	{
		return (this->location_ == cl_e::zi_trunk ||
			this->location_ == cl_e::zi_branch);
	}

	// --------------
	// constructor(s)
	// --------------
	Concealed(pgt::id_or_t fl_b, pgt::id_or_t cn_b, bounds_t bounds,
		  pt::idx_t fl_idx, cl_e loc, pt::idx_t loc_idx, route_e rt)
	    : VertexBase(pc::INVALID_IDX, vf_e::concealed), fl_idx_(fl_idx),
	      location_(loc), loc_idx_(loc_idx), fl_b_(fl_b), cn_b_(cn_b),
	      bounds_(bounds), route_(rt)
	{}

public:
	// --------------
	// factory fns
	// --------------
	static Concealed create(pgt::id_or_t fl_b, pgt::id_or_t cn_b,
				bounds_t bounds, pt::idx_t fl_idx, cl_e loc,
				pt::idx_t loc_idx, route_e rt)
	{
		return Concealed(fl_b, cn_b, bounds, fl_idx, loc, loc_idx, rt);
	}

	static Concealed parse(const route_params_t &rp)
	{
		pt::idx_t fl_idx{pc::INVALID_IDX};
		pvst::cl_e loc{cl_e::undefined};
		pt::idx_t loc_idx{pc::INVALID_IDX};
		bounds_t bounds{INVALID_BOUNDS};

		auto [fl_b, cn_b, r] = rp;
		return Concealed(fl_b, cn_b, bounds, fl_idx, loc, loc_idx, r);
	}

	// -------
	// getters
	// -------
	pt::idx_t get_fl_idx() const
	{
		return this->fl_idx_;
	}

	cl_e get_sl_type() const
	{
		return this->location_;
	}

	pt::idx_t get_sl_st_idx() const
	{
		return this->loc_idx_;
	}

	pgt::id_or_t get_fl_b() const
	{
		return this->fl_b_;
	}

	pgt::id_or_t get_cn_b() const
	{
		return this->cn_b_;
	}

	bounds_t get_bounds() const
	{
		return this->bounds_;
	}

	std::optional<route_params_t> get_route_params() const override
	{
		return route_params_t{this->get_fl_b(), this->get_cn_b(),
				      this->route_};
	}

	// ------
	// others
	// ------
	std::string as_str() const override
	{
		if (with_ai()) { // formed with a
			return pv_cmp::format("{}{}", this->fl_b_.as_str(),
					      this->cn_b_.as_str());
		}
		else { // formed with z
			return pv_cmp::format("{}{}", this->cn_b_.as_str(),
					      this->fl_b_.as_str());
		}
	}
};

class Smothered : public VertexBase
{
	pt::idx_t cn_idx_;    // idx of the concealed vertex
	pt::idx_t sm_st_idx_; // idx in the spanning tree for smothered vertex

	// b for boundary
	pgt::id_or_t cn_b_; // g or s
	pgt::id_or_t sm_b_; // e or w

	// is true when cn_b_ is an ancestor of sm_b_
	bool cn_b_is_ans_;

	cb_e cn_type_; // type of the concealed vertex (g or s)

	bounds_t bounds_;
	route_e route_;

	// --------------
	// constructor(s)
	// --------------
	Smothered(pgt::id_or_t cn_b, pgt::id_or_t sm_b, pt::idx_t cn_idx,
		  bool cn_b_is_ans, pt::idx_t sm_st_idx, cb_e sm_type,
		  bounds_t bounds, route_e route)
	    : VertexBase(pc::INVALID_IDX, vf_e::smothered), cn_idx_(cn_idx),
	      sm_st_idx_(sm_st_idx), cn_b_(cn_b), sm_b_(sm_b),
	      cn_b_is_ans_(cn_b_is_ans), cn_type_(sm_type), bounds_(bounds),
	      route_(route)
	{}

public:
	// -----------
	// factory fns
	// -----------
	static Smothered create(pgt::id_or_t cn_b, pgt::id_or_t sm_b,
				pt::idx_t cn_idx, bool cn_b_is_ans,
				pt::idx_t sm_st_idx, cb_e sm_type,
				bounds_t bounds, route_e route)
	{
		return Smothered(cn_b, sm_b, cn_idx, cn_b_is_ans, sm_st_idx,
				 sm_type, bounds, route);
	}

	static Smothered parse(const route_params_t &rp)
	{
		pt::idx_t cn_idx{pc::INVALID_IDX};
		pt::idx_t sm_st_idx{pc::INVALID_IDX};
		bool cn_b_is_ans{false};
		cb_e cn_type{cb_e::g};
		bounds_t bounds{INVALID_BOUNDS};

		auto [sm_b, cn_b, r] = rp;
		return Smothered(cn_b, sm_b, cn_idx, cn_b_is_ans, sm_st_idx,
				 cn_type, bounds, r);
	}

	// -------
	// getters
	// -------
	pt::idx_t get_cn_idx() const
	{
		return this->cn_idx_;
	}

	pt::idx_t get_sm_st_idx() const
	{
		return this->sm_st_idx_;
	}

	bounds_t get_bounds() const
	{
		return this->bounds_;
	}

	cb_e get_cn_type() const
	{
		return this->cn_type_;
	}

	bool is_cn_b_ancestor() const
	{
		return this->cn_b_is_ans_;
	}

	std::optional<route_params_t> get_route_params() const override
	{
		return route_params_t{this->sm_b_, this->cn_b_, this->route_};
	}

	std::string as_str() const override
	{
		if (this->cn_type_ == cb_e::g) {  // g
			if (this->cn_b_is_ans_) { // cn_b is ancestor of sm_b
				return pv_cmp::format("{}{}",
						      this->cn_b_.as_str(),
						      this->sm_b_.as_str());
			}
			else { // sm_b is ancestor of cn_b
				return pv_cmp::format("{}{}",
						      this->sm_b_.as_str(),
						      this->cn_b_.as_str());
			}
		}
		else { // s
			return pv_cmp::format("{}{}", this->sm_b_.as_str(),
					      this->cn_b_.as_str());
		}
	}
};

class MidiBubble : public VertexBase
{
	pt::idx_t g_cn_idx_;
	pt::idx_t s_cn_idx_;
	pgt::id_or_t g_;
	pgt::id_or_t s_;
	route_e route_;

	// --------------
	// constructor(s)
	// --------------
	MidiBubble(pt::idx_t g_cn_idx, pgt::id_or_t g, pt::idx_t s_cn_idx,
		   pgt::id_or_t s, route_e route)
	    : VertexBase(pc::INVALID_IDX, vf_e::midi), g_cn_idx_(g_cn_idx),
	      s_cn_idx_(s_cn_idx), g_(g), s_(s), route_(route)
	{}

public:
	// -----------
	// factory fns
	// -----------
	static MidiBubble create(pt::idx_t g_cn_idx, pgt::id_or_t g,
				 pt::idx_t s_cn_idx, pgt::id_or_t s,
				 route_e route)
	{
		// Validate indices at creation time
		if (g_cn_idx == pc::INVALID_IDX ||
		    s_cn_idx == pc::INVALID_IDX) {
			throw std::runtime_error("Cannot create MidiBubble "
						 "with invalid indices");
		}
		return MidiBubble(g_cn_idx, g, s_cn_idx, s, route);
	}

	static MidiBubble parse(const route_params_t &rp)
	{
		pt::idx_t g_cn_idx{pc::INVALID_IDX};
		pt::idx_t s_cn_idx{pc::INVALID_IDX};
		auto [g, s, r] = rp;
		return MidiBubble(g_cn_idx, g, s_cn_idx, s, r);
	}

	// -------
	// getters
	// -------
	bounds_t get_bounds() const
	{
		// Safety check - should never happen
		if (g_cn_idx_ == pc::INVALID_IDX ||
		    s_cn_idx_ == pc::INVALID_IDX) {
			return INVALID_BOUNDS;
		}
		return bounds_t{std::min(g_cn_idx_, s_cn_idx_),
				std::max(g_cn_idx_, s_cn_idx_)};
	}

	pgt::id_or_t get_g() const
	{
		return this->g_;
	}

	pgt::id_or_t get_s() const
	{
		return this->s_;
	}

	std::optional<route_params_t> get_route_params() const override
	{
		return route_params_t{this->get_g(), this->get_s(),
				      this->route_};
	}

	std::string as_str() const override
	{
		return pv_cmp::format("{}{}", this->g_.as_str(),
				      this->s_.as_str());
	}
};

// always has a dummy root vertex
class Tree
{
	std::vector<std::unique_ptr<pvst::VertexBase>> vertices;
	std::vector<pt::idx_t> parent_v; // parent of each vertex
	std::vector<std::vector<pt::idx_t>>
		children_v;  // children of each vertex
	pt::idx_t root_idx_; // index of the root vertex in the vertices vector

public:
	// --------------
	// constructor(s)
	// --------------

	Tree() : root_idx_(pc::INVALID_IDX)
	{}

	Tree(pt::idx_t expected_size) : Tree()
	{
		vertices.reserve(expected_size);
		this->parent_v = std::vector<pt::idx_t>(expected_size + 1,
							pc::INVALID_ID);
		this->children_v =
			std::vector<std::vector<pt::idx_t>>(1 + expected_size);
	}

	// ---------
	// getter(s)
	// ---------

	pt::idx_t vtx_count() const
	{
		return this->vertices.size();
	}

	pt::idx_t root_idx() const
	{
		return this->root_idx_;
	}

	const pvst::VertexBase &get_root() const
	{
		return this->get_vertex(this->root_idx());
	}

	const pvst::VertexBase &get_vertex(pt::idx_t v_idx) const
	{
		return *this->vertices[v_idx];
	}

	const pvst::VertexBase *get_vertex_const_ptr(pt::idx_t v_idx) const
	{
		return this->vertices[v_idx].get();
	}

	pvst::VertexBase &get_vertex_mut(pt::idx_t v_idx)
	{
		return *this->vertices[v_idx];
	}

	const pvst::VertexBase &get_parent(pt::idx_t v_idx) const
	{
		return *this->vertices[parent_v[v_idx]];
	}

	pt::idx_t get_parent_idx(pt::idx_t v_idx) const
	{
		return this->parent_v[v_idx];
	}

	// TODO: rename to get_children_idxs()
	const std::vector<pt::idx_t> &get_children(pt::idx_t v_idx) const
	{
		return this->children_v[v_idx];
	}

	bool is_leaf(pt::idx_t v_idx) const
	{
		return v_idx >= this->children_v.size() ||
		       this->children_v[v_idx].empty();
	}

	// ---------
	// setter(s)
	// ---------

	/**
	 * @brief Compute the heights of all vertices in the tree
	 */
	void comp_heights()
	{
		if (this->root_idx_ == pc::INVALID_IDX) {
			throw std::logic_error("Root index is not set");
		}

		// reset heights
		for (auto &v : this->vertices) {
			v->set_height(0);
		}

		// compute heights using DFS
		std::stack<pt::idx_t> s;
		s.push(this->root_idx_);

		while (!s.empty()) {
			pt::idx_t v_idx = s.top();
			s.pop();

			for (pt::idx_t child_v_idx :
			     this->get_children(v_idx)) {
				auto &child_v =
					this->get_vertex_mut(child_v_idx);
				child_v.set_height(
					this->get_vertex(v_idx).get_height() +
					1);
				s.push(child_v_idx);
			}
		}
	}

	void set_root_idx(pt::idx_t v_idx)
	{
		if (this->root_idx_ != pc::INVALID_IDX) {
			throw std::logic_error("Root index is already set");
		}

		if (v_idx >= this->vertices.size()) {
			throw std::out_of_range("Vertex index out of range");
		}

		this->root_idx_ = v_idx;
	}

	/**
	 * @brief Add a vertex to the tree
	 * @param v Vertex to be added
	 * @return Index of the added vertex
	 */
	template <typename T> pt::idx_t add_vertex(T v)
	{
		pt::idx_t v_idx = this->vertices.size();
		v.set_idx(v_idx); // set the index of the vertex

		while (v_idx >= this->parent_v.size()) {
			this->parent_v.push_back(
				pc::INVALID_ID); // ensure parent_v has enough
						 // space
		}

		while (v_idx >= this->children_v.size()) {
			this->children_v.push_back(
				std::vector<pt::idx_t>()); // ensure children_v
							   // has enough space
		}

		// create a unique pointer to the vertex and add it to the
		// vertices vector
		auto ptr = std::make_unique<T>(v);
		this->vertices.push_back(std::move(ptr));

		return v_idx;
	}

	void add_edge(pt::idx_t parent, pt::idx_t child)
	{
		while (child >= this->parent_v.size()) {
			this->parent_v.push_back(pc::INVALID_ID);
		}
		this->parent_v[child] = parent;
		while (parent >= this->children_v.size()) {
			this->children_v.push_back(std::vector<pt::idx_t>());
		}
		this->children_v[parent].push_back(child);
	}

	void del_edge(pt::idx_t parent, pt::idx_t child)
	{
		this->parent_v[child] = pc::INVALID_IDX;

		std::vector<pt::idx_t> &children = this->children_v[parent];
		auto it = std::find(children.begin(), children.end(), child);
		if (it != children.end()) {
			children.erase(it);
		}
	}

	// ----
	// misc
	// ----
	void print_dot() const
	{
		std::cout << pv_cmp::format("graph G {{\n"
					    "\trankdir = TD;\n"
					    "\tnode[shape = circle];\n"
					    "\tedge [arrowhead=vee];\n");

		// print vertices
		for (pt::idx_t i{}; i < this->vtx_count(); i++) {
			std::cout
				<< pv_cmp::format("\t{} [label=\"{}\"];\n", i,
						  this->get_vertex(i).as_str());
		}

		// print edges
		for (pt::idx_t i{}; i < this->vtx_count(); i++) {
			for (pt::idx_t c : this->get_children(i)) {
				std::cout << pv_cmp::format("\t{} -- {};\n", i,
							    c);
			}
		}

		std::cout << "}" << std::endl;
	}
};

} // namespace povu::pvst

#endif

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pvst = povu::pvst;
