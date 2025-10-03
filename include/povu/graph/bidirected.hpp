#ifndef BIDIRECTED_HPP
#define BIDIRECTED_HPP

#include <algorithm>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <iterator>
#include <map>
#include <optional>
#include <ostream>
#include <queue>
#include <set>
#include <stack>
#include <string>
#include <string_view>
#include <sys/types.h>
#include <tuple>
#include <unordered_set>
#include <vector>

#include <liteseq/gfa.h>

"#include "povu/common/compat.hpp"
"#include "povu/common/constants.hpp"
"#include "povu/common/log.hpp"
"#include "povu/common/utils.hpp"
#include "refs/refs.hpp"
#include "types.hpp"

namespace povu::bidirected
{
inline constexpr std::string_view MODULE = "povu::bidirected";

using namespace povu::types::graph;
namespace pgt = povu::types::graph;

namespace lq = liteseq;

// undirected edge
// stores the index of the vertex in the graph not the id
class Edge
{
	pt::idx_t v1_idx_;
	pgt::v_end_e v1_end_;
	pt::idx_t v2_idx_;
	pgt::v_end_e v2_end_;

public:
	// --------------
	// constructor(s)
	// --------------
	Edge(pt::idx_t v1_idx, pgt::v_end_e v1_end, pt::idx_t v2_idx,
	     pgt::v_end_e v2_end);

	// ---------
	// getter(s)
	// ---------
	pt::idx_t get_v1_idx() const;
	pt::idx_t &get_v1_idx_mut();
	pgt::v_end_e get_v1_end() const;
	pt::idx_t get_v2_idx() const;
	pt::idx_t &get_v2_idx_mut();
	pgt::v_end_e get_v2_end() const;
	pgt::side_n_idx_t get_other_vtx(pt::idx_t v_id,
					pgt::v_end_e v_end) const;
	pgt::side_n_idx_t
	get_other_vtx(pt::idx_t v_id) const; // if you don't care for self loops
};

class Vertex
{
	pt::id_t v_id_;	    // this is the sequence name in the GFA file. Maybe
			    // Should support strings.
	std::string label_; // or sequence

	// indexes to the edge vector in Graph
	std::set<pt::idx_t> e_l;
	std::set<pt::idx_t> e_r;

public:
	// --------------
	// constructor(s)
	// --------------
	Vertex(pt::id_t v_id, const std::string &label = "");

	// ---------
	// getter(s)
	// ---------
	pt::id_t id() const;
	const std::string &get_label() const;
	std::string get_rc_label() const; // reverse complement of the label
	const std::set<pt::idx_t> &get_edges_l() const;
	const std::set<pt::idx_t> &get_edges_r() const;

	// ---------
	// setter(s)
	// ---------
	void add_edge_l(pt::idx_t e_idx);
	void add_edge_r(pt::idx_t e_idx);
};

class VariationGraph
{
	std::vector<Vertex> vertices;
	std::set<pgt::side_n_id_t> tips_; // the set of side and id of the tips
	pu::TwoWayMap<std::size_t, std::size_t>
		v_id_to_idx_; // TODO: reserve size

	std::vector<Edge> edges;

	// i is the vertex index, j is the ref index
	std::vector<std::vector<std::vector<pt::idx_t>>> vertex_to_step_matrix_;
	pr::Refs refs_;

	lq::gfa_props *gfa;

	// i is the vertex index, j is the ref index
	// std::vector<std::vector<std::vector<pt::idx_t>>>
	// vertex_to_step_matrix_; std::vector<pgt::ref_walk_t> ref_matrix_;
	// pr::Refs refs_ = pr::Refs(0); // has no refs by default

public:
	// --------------
	// constructor(s)
	// --------------
	VariationGraph(pt::idx_t vtx_count, pt::idx_t edge_count,
		       pt::idx_t ref_count);
	VariationGraph(lq::gfa_props *gfa);

	~VariationGraph()
	{
		if (this->gfa != nullptr) {
			gfa_free(this->gfa);
		}
	}

	// -----------------
	// factory method(s)
	// -----------------
	// return a vector of connected components as VG objects
	static std::vector<VariationGraph *>
	componetize(const VariationGraph &g);

	// ---------
	// getter(s)
	// ---------
	pt::idx_t v_id_to_idx(pt::id_t v_id) const;
	pt::id_t v_idx_to_id(pt::idx_t v_idx) const;

	pt::idx_t vtx_count() const;
	pt::idx_t edge_count() const;
	const std::set<pgt::side_n_id_t> &tips() const;
	const Edge &get_edge(pt::idx_t e_idx) const;
	Edge &get_edge_mut(pt::idx_t e_idx);
	// TODO replace vertex with v?
	const Vertex &get_vertex_by_idx(pt::idx_t v_idx) const;
	const Vertex &get_vertex_by_id(pt::id_t v_id) const;
	Vertex &get_vertex_mut_by_id(pt::id_t v_id);

	// ref
	const std::string &get_sample_name(pt::id_t ref_id) const;
	const pr::Ref &get_ref_by_id(pt::id_t ref_id) const;
	pr::Ref &get_ref_by_id_mut(pt::id_t ref_id);
	std::optional<pt::id_t> get_ref_id(std::string_view ref_tag) const;
	pt::id_t ref_count() const;
	bool has_refs() const;
	/**
	 * if PanSN the prefix will be in sample
	 */
	std::set<pt::id_t> get_shared_samples(pt::id_t ref_id) const;
	// sometimes the sample name is referred to as a prefix
	std::set<pt::id_t>
	get_refs_in_sample(std::string_view sample_name) const;
	const lq::ref *get_ref_vec(pt::id_t ref_id) const;
	pt::idx_t get_ref_count() const;
	const std::vector<pt::idx_t> &
	get_vertex_ref_idxs(pt::idx_t v_idx, pt::id_t ref_id) const;
	const std::vector<std::string> &get_genotype_col_names() const;
	std::vector<std::vector<std::string>> get_blank_genotype_cols() const;
	const pt::op_t<pt::idx_t> &get_ref_gt_col_idx(pt::id_t ref_id) const;

	// ---------
	// setter(s)
	// ---------
	void add_tip(pt::id_t v_id, pgt::v_end_e end);
	// returns the index (v_idx) of the added vertex
	pt::idx_t add_vertex(pt::id_t v_id, const std::string &label);
	// returns the index (e_idx) of the added edge
	pt::idx_t add_edge(pt::id_t v1_id, pgt::v_end_e v1_end, pt::id_t v2_id,
			   pgt::v_end_e v2_end);
	void add_all_refs(lq::ref **refs, pt::idx_t ref_count);
	// pt::id_t add_ref(const std::string &label, char delim);
	void shrink_to_fit();
	void set_vtx_ref_idx(pt::id_t v_id, pt::id_t ref_id,
			     pt::idx_t step_idx);
	void gen_genotype_metadata();

	// -----
	// other
	// -----
	void summary(bool print_tips) const;
	void print_dot(std::ostream &os) const;
};

typedef VariationGraph VG;
} // namespace povu::bidirected

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace bd = povu::bidirected;

#endif
