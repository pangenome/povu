#ifndef IT_SL_TREE_HPP
#define IT_SL_TREE_HPP

#include <cstdlib>
#include <map>
#include <optional>
#include <set>

#include <log.h>
#include <quilt/constants.hpp> // for
#include <quilt/shim.hpp>      // for format
#include <quilt/types.hpp>     // for qt

namespace ita::slice_tree
{

/**/
enum class comp_type : qt::u8 {
	NO_OVERLAP = 0,
	EXISTS = 1,
	REPLACE_ALT = 2,
	MERGE_EXTEND = 3,  // merge extend or merge left
	MERGE_REPLACE = 4, // merge and replace or merge right
	INSERT_ALT = 5,
	CONTAINED = 6,	// new is contained within what exists
	EXTEND_ALT = 7, // for inverted repeats, add extra alt
};

std::string to_string(comp_type ct);

// compare result options
// restult of comparing two intervals
struct comp_res {
	comp_type ct;
	// qt::u32 shift_ref{0};
	qt::u32 ai{pc::INVALID_IDX}; // alt index

	// ------------
	// constructors
	// ------------

	comp_res() = delete;

	comp_res(comp_type ct_) : ct{ct_}
	{}

	comp_res(comp_type ct_, qt::u32 alt_idx) : ct{ct_}, ai{alt_idx}
	{}
};

/* insert specific */
enum class child_type : qt::u8 {
	LEFT = 0,
	RIGHT = 1
};

struct insert_opts {
	qt::u32 parent_idx; // index of the parent node
	child_type ct;	    // add new node as left or right child of parent
};

enum class update_type : qt::u8 {
	DO_NOTHING = 1,
	REPLACE_ALT = 2,
	MERGE_EXTEND = 3,  // also, merge extend
	MERGE_REPLACE = 4, // merge and replace
	INSERT_LEAF = 5,   // always inserts a leaf node
	INSERT_ALT = 6,	   // adds an alt where none existed
	EXTEND_ALT = 7,	   // for invereted repeats, adds extra alts
};

std::string to_string(update_type ut);

struct update_params {
	update_type ut;
	qt::u32 existing_idx; // index of existing vertex in the tree to update
	qt::u32 alt_idx{pc::INVALID_IDX}; // index of the alt to update
	std::optional<insert_opts> ip;

	// ------------
	// constructors
	// ------------
	update_params() = delete;

	update_params(update_type ut_) : ut{ut_}, ip{std::nullopt}
	{}

	update_params(update_type ut_, insert_opts ip_) : ut{ut_}, ip{ip_}
	{}

	update_params(update_type ut_, qt::u32 ei)
	    : ut{ut_}, existing_idx{ei}, ip{std::nullopt}
	{}

	update_params(update_type ut_, qt::u32 ei, qt::u32 a_idx)
	    : ut{ut_}, existing_idx{ei}, alt_idx{a_idx}, ip{std::nullopt}
	{}
};

struct alt {
	qt::u32 h_idx;
	qt::u32 h_start;
	qt::u32 len;
};

struct vertex {
private:
	// --------------------
	// vertex specific data
	// --------------------

	qt::u32 r_h_start_;			  // ref hap start
	std::map<qt::u32, std::vector<alt>> alts; // alt hap idx -> alts
	// alt hap length -> haps that match the ref at this length
	std::map<qt::u32, std::set<qt::u32>> len_haps;

	// ------------------
	// interval tree data
	// ------------------

	qt::u32 parent_;
	qt::u32 left_;
	qt::u32 right_;

	// ----------------
	// internal methods
	// ----------------
	vertex() = default;

public:
	// static factory methods
	static vertex create_root(qt::u32 ref_h_start, qt::u32 alt_h_idx,
				  qt::u32 alt_h_start, qt::u32 len)
	{
		vertex v;
		v.r_h_start_ = ref_h_start;
		v.parent_ = pc::INVALID_IDX;
		v.left_ = pc::INVALID_IDX;
		v.right_ = pc::INVALID_IDX;

		v.insert_alt(alt_h_idx, alt_h_start, len);

		return v;
	}

	static vertex create_leaf(qt::u32 ref_h_start, qt::u32 alt_h_idx,
				  qt::u32 alt_h_start, qt::u32 len,
				  qt::u32 parent_idx)
	{
		vertex v;
		v.r_h_start_ = ref_h_start;
		v.parent_ = parent_idx;
		v.left_ = pc::INVALID_IDX;
		v.right_ = pc::INVALID_IDX;

		v.insert_alt(alt_h_idx, alt_h_start, len);

		return v;
	}

	// -------
	// getters
	// -------

	[[nodiscard]]
	qt::u32 get_r_start() const
	{
		return this->r_h_start_;
	}

	[[nodiscard]]
	const std::map<qt::u32, std::set<qt::u32>> get_same_len_alts() const
	{
		return this->len_haps;
	}

	[[nodiscard]]
	const std::set<qt::u32> *get_len_alts(qt::u32 len) const
	{
		if (qs::contains(this->len_haps, len))
			return &this->len_haps.at(len);
		else
			return nullptr;
	}

	[[nodiscard]]
	qt::u32 left() const
	{
		return this->left_;
	}

	[[nodiscard]]
	qt::u32 right() const
	{
		return this->left_;
	}

	[[nodiscard]]
	qt::u32 parent() const
	{
		return this->parent_;
	}

	[[nodiscard]]
	bool only_right() const
	{
		return this->left_ == pc::INVALID_IDX &&
		       this->right_ != pc::INVALID_IDX;
	}

	[[nodiscard]]
	bool only_left() const
	{
		return this->left_ != pc::INVALID_IDX &&
		       this->right_ == pc::INVALID_IDX;
	}

	[[nodiscard]]
	bool is_leaf() const
	{
		return this->left_ == pc::INVALID_IDX &&
		       this->right_ == pc::INVALID_IDX;
	}

	[[nodiscard]]
	bool has_alts() const
	{
		return !this->alts.empty();
	}

	[[nodiscard]]
	const std::vector<alt> &get_alts(qt::u32 alt_h_idx) const
	{
		return this->alts.at(alt_h_idx);
	}

	// ---------
	// modifiers
	// ---------

	[[nodiscard]]
	std::vector<alt> &get_alts_mut(qt::u32 alt_h_idx)
	{
		return this->alts.at(alt_h_idx);
	}

	void set_parent(qt::u32 p_idx)
	{
		if (p_idx == this->get_r_start()) {
			log_fatal("Setting parent to self");
			std::exit(EXIT_FAILURE);
		}

		this->parent_ = p_idx;
	}

	void set_left(qt::u32 l_idx)
	{
		if (this->left() != pc::INVALID_IDX) {
			log_fatal("Adding leaf to non-empty left child");
			std::exit(EXIT_FAILURE);
		}

		this->left_ = l_idx;
	}

	void set_right(qt::u32 r_idx)
	{
		if (this->right() != pc::INVALID_IDX) {
			log_fatal("Adding leaf to non-empty right child");
			std::exit(EXIT_FAILURE);
		}

		this->right_ = r_idx;
	}

	[[nodiscard]]
	bool alts_has_hap(qt::u32 alt_h_idx) const
	{
		return qs::contains(this->alts, alt_h_idx);
	}

	[[nodiscard]]
	comp_res check_update_type(qt::u32 ref_h_start, qt::u32 alt_h_idx,
				   qt::u32 alt_h_start, qt::u32 len)
	{

		auto check_alt = [&](const alt &a) -> std::optional<comp_type>
		{
			if (a.h_start == alt_h_start && a.len == len)
				return comp_type::EXISTS; // already exists

			/* containment */

			// new is contained within what exists
			if (a.h_start < alt_h_start &&
			    a.h_start + a.len >= alt_h_start + len) {
				return comp_type::CONTAINED;
			}

			// what exists is contained within new
			if (alt_h_start < a.h_start &&
			    alt_h_start + len >= a.h_start + a.len) {
				return comp_type::REPLACE_ALT;
			}

			/* overlapping and need to merge*/

			if (a.h_start < alt_h_start &&
			    a.h_start + a.len >= alt_h_start &&
			    alt_h_start + len >= a.h_start + a.len) {
				return comp_type::MERGE_EXTEND;
			}

			if (alt_h_start < a.h_start &&
			    alt_h_start + len >= a.h_start &&
			    a.h_start + a.len >= alt_h_start + len) {
				return comp_type::MERGE_REPLACE;
			}

			return std::nullopt;
		};

		auto check_alts = [&](comp_type fallback) -> comp_res
		{
			qt::u32 N = this->get_alts(alt_h_idx).size();
			for (qt::u32 i{}; i < N; i++) {
				const alt &a = this->get_alts(alt_h_idx).at(i);
				std::optional<comp_type> opt_ct = check_alt(a);

				if (opt_ct.has_value())
					return {opt_ct.value(), i};
			}

			return {fallback};
		};

		if (this->get_r_start() == ref_h_start &&
		    !this->alts_has_hap(alt_h_idx))
			return {comp_type::INSERT_ALT};

		if (this->get_r_start() == ref_h_start)
			return check_alts(comp_type::EXTEND_ALT);
		else // when the ref starts are not equal
			return check_alts(comp_type::NO_OVERLAP);
	}

	// used when inserting a new alt hap that does not exist yet
	void insert_alt(qt::u32 alt_h_idx, qt::u32 alt_h_start, qt::u32 len)
	{
		alt a;
		a.h_idx = alt_h_idx;
		a.h_start = alt_h_start;
		a.len = len;

		this->alts.insert({alt_h_idx, {a}});
		this->len_haps[len].insert(alt_h_idx);
	}

	// used when adding extra alts for inverted repeats
	void append_alt(qt::u32 alt_h_idx, qt::u32 alt_h_start, qt::u32 len)
	{
		alt a;
		a.h_idx = alt_h_idx;
		a.h_start = alt_h_start;
		a.len = len;

		this->alts[alt_h_idx].emplace_back(a);
		this->len_haps[len].insert(alt_h_idx);
	}

	void remove_alt(qt::u32 alt_h_idx, qt::u32 alt_idx)
	{
		std::vector<alt> &alts_vec = this->get_alts_mut(alt_h_idx);
		if (alt_idx >= alts_vec.size()) {
			std::string err = qs::format(
				"Alt index {} out of bounds for alt hap {}",
				alt_idx, alt_h_idx);
			log_fatal("%s", err.c_str());
			std::exit(EXIT_FAILURE);
		}

		qt::u32 len = alts_vec.at(alt_idx).len;
		alts_vec.erase(alts_vec.begin() + alt_idx);

		this->len_haps[len].erase(alt_h_idx);
		if (this->len_haps[len].empty())
			this->len_haps.erase(len);
	}

	void replace_alt(qt::u32 alt_h_idx, qt::u32 alt_h_start, qt::u32 len,
			 qt::u32 a_idx)
	{
		this->remove_alt(alt_h_idx, a_idx); // remove existing alt
		this->insert_alt(alt_h_idx, alt_h_start, len); // add new alt
	}

	void extend_alt(qt::u32 alt_h_idx, qt::u32 alt_h_start, qt::u32 len,
			qt::u32 a_idx)
	{
		alt &a = this->get_alts_mut(alt_h_idx).at(a_idx);
		qt::u32 old_len = a.len;
		qt::u32 new_len = (alt_h_start + len) - a.h_start;

		a.h_start = alt_h_start;
		a.len = new_len;

		// update len_haps map
		this->len_haps[old_len].erase(alt_h_idx);
		this->len_haps[new_len].insert(alt_h_idx);
	}
};

// slice tree, range tree, interval tree, segment tree etc
struct slice_tree {
private:
	qt::u32 ref_h_idx_{pc::INVALID_IDX}; // ref hap index
	qt::u32 root_{pc::INVALID_IDX};
	std::map<qt::u32, vertex> vertices{};

	void set_root(qt::u32 r_idx)
	{
		this->root_ = r_idx;
	}

	void add_root(qt::u32 ref_h_start, qt::u32 alt_h_idx,
		      qt::u32 alt_h_start, qt::u32 len)
	{
		// std::cerr << "Adding root vertex at ref_h_start " <<
		// ref_h_start
		//	  << "\n";

		vertex v = vertex::create_root(ref_h_start, alt_h_idx,
					       alt_h_start, len);

		this->set_root(ref_h_start);
		this->vertices.emplace(ref_h_start, v);
	}

	void add_leaf(qt::u32 ref_h_start, qt::u32 alt_h_idx,
		      qt::u32 alt_h_start, qt::u32 len, const update_params &up)
	{
		vertex v =
			vertex::create_leaf(ref_h_start, alt_h_idx, alt_h_start,
					    len, up.ip->parent_idx);

		this->vertices.emplace(ref_h_start, v);

		/* update the parent */
		vertex &parent = this->get_vertex_mut(up.ip->parent_idx);
		if (up.ip->ct == child_type::LEFT)
			parent.set_left(ref_h_start);
		else // RIGHT
			parent.set_right(ref_h_start);
	}

	void insert_alt(qt::u32 alt_h_idx, qt::u32 alt_h_start, qt::u32 len,
			const update_params &up)
	{
		vertex &v = this->get_vertex_mut(up.existing_idx);
		v.insert_alt(alt_h_idx, alt_h_start, len);
	}

	void append_alt(qt::u32 alt_h_idx, qt::u32 alt_h_start, qt::u32 len,
			const update_params &up)
	{
		vertex &v = this->get_vertex_mut(up.existing_idx);
		v.append_alt(alt_h_idx, alt_h_start, len);
	}

	void replace_alt(qt::u32 alt_h_idx, qt::u32 alt_h_start, qt::u32 len,
			 const update_params &up)
	{
		vertex &v = this->get_vertex_mut(up.existing_idx);
		v.replace_alt(alt_h_idx, alt_h_start, len, up.alt_idx);
	}

	void merge_extend(qt::u32 alt_h_idx, qt::u32 alt_h_start, qt::u32 len,
			  const update_params &up)
	{
		vertex &v = this->get_vertex_mut(up.existing_idx);
		v.extend_alt(alt_h_idx, alt_h_start, len, up.alt_idx);
	}

	void remove_leaf(qt::u32 v_idx)
	{
		qt::u32 root_v_idx = this->root();
		if (v_idx == root_v_idx) {
			// removing root
			this->set_root(pc::INVALID_IDX);
			this->vertices.erase(v_idx);
			return;
		}

		// update parent
		vertex &v = this->get_vertex_mut(v_idx);
		qt::u32 p_idx = v.parent();
		if (p_idx == pc::INVALID_IDX) {
			std::string err = qs::format(
				"Non root leaf vertex {} has no parent", v_idx);
			log_fatal("%s", err.c_str());
			std::exit(EXIT_FAILURE);
		}

		vertex &parent = this->get_vertex_mut(p_idx);
		if (parent.left() == v_idx)
			parent.set_left(pc::INVALID_IDX);
		else if (parent.right() == v_idx)
			parent.set_right(pc::INVALID_IDX);

		// remove vertex
		this->vertices.erase(v_idx);
	}

	[[nodiscard]]
	qt::u32 left_most_leaf(qt::u32 left_most)
	{
		while (left_most != pc::INVALID_IDX) {
			qt::u32 lc = this->get_vertex_mut(left_most).left();
			if (lc == pc::INVALID_IDX) {
				break;
			}

			left_most = lc;
		}

		return left_most;
	}

	void remove_vertex(qt::u32 v_idx)
	{
		vertex &v = this->get_vertex_mut(v_idx);

		if (v.is_leaf())
			return this->remove_leaf(v_idx);

		qt::u32 p_idx = v.parent();

		if (v.only_left()) {
			// update the left child
			qt::u32 left_idx = v.left();
			this->get_vertex_mut(left_idx).set_parent(p_idx);

			// update parent
			if (p_idx != pc::INVALID_IDX) {
				vertex &parent = this->get_vertex_mut(p_idx);
				parent.set_left(left_idx);
			}
			else {
				this->set_root(left_idx); // parent was root
			}

			// remove vertex
			this->vertices.erase(v_idx);
			return;
		}

		/* has a right child, update left_most to take v's place */

		qt::u32 m = this->left_most_leaf(v.right());
		vertex &left_most = this->get_vertex_mut(m);

		// update left_most's exising parent
		this->get_vertex_mut(left_most.parent())
			.set_left(pc::INVALID_IDX);

		// move left_most to v's position
		left_most.set_parent(v.parent());

		if (v.left() != pc::INVALID_IDX) {
			this->get_vertex_mut(v.left()).set_parent(m);
			left_most.set_left(v.left());
		}

		if (v.right() != pc::INVALID_IDX) {
			this->get_vertex_mut(v.right()).set_parent(m);
			left_most.set_right(v.right());
		}
	}

	void merge_replace(qt::u32 ref_h_start, qt::u32 alt_h_idx,
			   qt::u32 alt_h_start, const update_params &up)
	{
		vertex &v = this->get_vertex_mut(up.existing_idx);

		qt::u32 a_idx = up.alt_idx;

		const std::vector<alt> &alts_vec = v.get_alts(alt_h_idx);
		alt a = alts_vec.at(a_idx);

		v.remove_alt(alt_h_idx, a_idx);

		// if the vertex ends up with no alts, remove the vertex
		if (!v.has_alts())
			this->remove_vertex(up.existing_idx);

		qt::u32 new_len = (a.h_start + a.len) - alt_h_start;

		this->add_vertex(ref_h_start, alt_h_idx, alt_h_start, new_len);
	}

	update_params find_insert_point_traverse(qt::u32 ref_h_start,
						 qt::u32 alt_h_idx,
						 qt::u32 alt_h_start,
						 qt::u32 len)
	{
		// INFO("Finding insert point for ref_h_start {} size {}",
		//      ref_h_start, this->size());

		std::set<qt::u32> visited;
		qt::u32 ctr{};

		qt::u32 curr_v_idx = this->root_;
		while (true) {

			// std::cerr << "Ctr " << ctr << " v " << curr_v_idx
			//	  << "\n";

			if (ctr++ >= this->size()) {
				log_fatal("Exceeded max traversal steps");
				std::exit(EXIT_FAILURE);
			}

			if (qs::contains(visited, curr_v_idx)) {
				log_fatal("Revist vertex %ul", curr_v_idx);
				std::exit(EXIT_FAILURE);
			}

			visited.insert(curr_v_idx);

			vertex &curr = this->get_vertex_mut(curr_v_idx);
			if (curr.alts_has_hap(alt_h_idx)) {
				comp_res cr = curr.check_update_type(
					ref_h_start, alt_h_idx, alt_h_start,
					len);

				switch (cr.ct) {
				case comp_type::NO_OVERLAP:
					break; // continue tree traversal
				case comp_type::CONTAINED:
				case comp_type::EXISTS:
					return {update_type::DO_NOTHING};
				case comp_type::INSERT_ALT:
					return {update_type::INSERT_ALT,
						curr_v_idx};
				case comp_type::REPLACE_ALT:
					return {update_type::REPLACE_ALT,
						curr_v_idx, cr.ai};
				case comp_type::MERGE_EXTEND:
					return {update_type::MERGE_EXTEND,
						curr_v_idx, cr.ai};
				case comp_type::MERGE_REPLACE:
					return {update_type::MERGE_REPLACE,
						curr_v_idx, cr.ai};
				case comp_type::EXTEND_ALT:
					log_fatal("Did not expect EXTEND_ALT "
						  "here");
					std::exit(EXIT_FAILURE);
				}
			}

			if (ref_h_start < curr.get_r_start()) {
				if (curr.left() == pc::INVALID_IDX) {
					return {update_type::INSERT_LEAF,
						{curr_v_idx, child_type::LEFT}};
				}

				curr_v_idx = curr.left();
			}
			else if (ref_h_start > curr.get_r_start()) {
				if (curr.left() == pc::INVALID_IDX) {
					return {update_type::INSERT_LEAF,
						{curr_v_idx,
						 child_type::RIGHT}};
				}

				curr_v_idx = curr.right();
			}

			if (ref_h_start == curr.get_r_start()) {
				// should not reach here
				log_fatal("Reached unexpected code path");
				std::exit(EXIT_FAILURE);
			}
		}
	}

	// for existing start in ref_h_start
	update_params handle_existing_start(qt::u32 ref_h_start,
					    qt::u32 alt_h_idx,
					    qt::u32 alt_h_start, qt::u32 len)
	{
		vertex &curr = this->get_vertex_mut(ref_h_start);
		comp_res cr = curr.check_update_type(ref_h_start, alt_h_idx,
						     alt_h_start, len);

		switch (cr.ct) {
		case comp_type::CONTAINED:
		case comp_type::EXISTS:
			return {update_type::DO_NOTHING};
		case comp_type::REPLACE_ALT:
			return {update_type::REPLACE_ALT, ref_h_start, cr.ai};
		case comp_type::MERGE_EXTEND:
			return {update_type::MERGE_EXTEND, ref_h_start, cr.ai};
		case comp_type::MERGE_REPLACE:
			return {update_type::MERGE_REPLACE, ref_h_start, cr.ai};
		case comp_type::INSERT_ALT:
			return {update_type::INSERT_ALT, ref_h_start};
		case comp_type::EXTEND_ALT:
			return {update_type::EXTEND_ALT, ref_h_start};
		case comp_type::NO_OVERLAP:
			std::string err =
				qs::format("Did not expect comp type {}",
					   to_string(cr.ct));
			log_fatal("%s", err.c_str());
			std::exit(EXIT_FAILURE);
		}

		std::string err = qs::format("Reached unexpected code path {}",
					     to_string(cr.ct));
		log_fatal("%s", err.c_str());
		std::exit(EXIT_FAILURE);
	}

	update_params find_insert_point(qt::u32 ref_h_start, qt::u32 alt_h_idx,
					qt::u32 alt_h_start, qt::u32 len)
	{
		if (qs::contains(vertices, ref_h_start)) {
			// is contained somehow
			return handle_existing_start(ref_h_start, alt_h_idx,
						     alt_h_start, len);
		}

		return this->find_insert_point_traverse(ref_h_start, alt_h_idx,
							alt_h_start, len);
	}

	void add_internal(qt::u32 ref_h_start, qt::u32 alt_h_idx,
			  qt::u32 alt_h_start, qt::u32 len)
	{
		update_params up = this->find_insert_point(
			ref_h_start, alt_h_idx, alt_h_start, len);

		switch (up.ut) {
		case update_type::DO_NOTHING:
			return;
		case update_type::INSERT_LEAF:
			return this->add_leaf(ref_h_start, alt_h_idx,
					      alt_h_start, len, up);
		case update_type::INSERT_ALT:
			return this->insert_alt(alt_h_idx, alt_h_start, len,
						up);
		case update_type::EXTEND_ALT:
			return this->append_alt(alt_h_idx, alt_h_start, len,
						up);
		case update_type::REPLACE_ALT:
			return this->replace_alt(alt_h_idx, alt_h_start, len,
						 up);
		case update_type::MERGE_EXTEND:
			return this->merge_extend(alt_h_idx, alt_h_start, len,
						  up);
		case update_type::MERGE_REPLACE:
			return this->merge_replace(ref_h_start, alt_h_idx,
						   alt_h_start, up);
		}
	}

public:
	// --------------
	// constructor(s)
	// --------------
	slice_tree() = delete;

	slice_tree(qt::u32 ref_h_idx) : ref_h_idx_{ref_h_idx}
	{}

	// -------
	// getters
	// -------

	[[nodiscard]]
	qt::u32 get_ref_hap_idx() const
	{
		return this->ref_h_idx_;
	}

	[[nodiscard]]
	qt::u32 root() const
	{
		return this->root_;
	}

	[[nodiscard]]
	qt::u32 size() const
	{
		return this->vertices.size();
	}

	[[nodiscard]]
	bool is_empty() const
	{
		{
#ifdef DEBUG // Check for inconsistency
			if (this->vertices.empty() &&
			    this->root() != pc::INVALID_IDX)
				throw std::runtime_error(
					"Interval tree inconsistency: root "
					"index set "
					"but vertices map is empty");

			if (!this->vertices.empty() &&
			    this->root() == pc::INVALID_IDX)
				throw std::runtime_error(
					"Interval tree inconsistency: vertices "
					"map "
					"is not empty but root index is "
					"invalid");
#endif
		}
		return this->root() == pc::INVALID_IDX &&
		       this->vertices.empty();
	}

	[[nodiscard]]
	const std::map<qt::u32, vertex> &get_vertices() const
	{
		return this->vertices;
	}

	// ---------
	// modifiers
	// ---------

	vertex &get_vertex_mut(qt::u32 v_idx)
	{
		return this->vertices.at(v_idx);
	}

	void add_vertex(qt::u32 ref_h_start, qt::u32 alt_h_idx,
			qt::u32 alt_h_start, qt::u32 len)
	{
		if (this->is_empty())
			this->add_root(ref_h_start, alt_h_idx, alt_h_start,
				       len);
		else
			this->add_internal(ref_h_start, alt_h_idx, alt_h_start,
					   len);
	}
};

using st = slice_tree;

} // namespace ita::slice_tree

namespace ist = ita::slice_tree;

#endif // IT_SL_TREE_HPP
