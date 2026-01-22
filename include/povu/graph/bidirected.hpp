#ifndef BIDIRECTED_HPP
#define BIDIRECTED_HPP

#include <cstddef>	 // for size_t
#include <iostream>	 // for ostream
#include <liteseq/gfa.h> // for gfa_free, gfa_props
#include <optional>	 // for optional
#include <set>		 // for set
#include <string>	 // for string, basic_string
#include <string_view>	 // for string_view
#include <vector>	 // for vector

#include "liteseq/refs.h"	 // for ref
#include "povu/common/core.hpp"	 // for pt, idx_t, id_t, op_t
#include "povu/common/utils.hpp" // for pu, TwoWayMap
#include "povu/graph/types.hpp"	 // for v_end_e, side_n_id_t, side_n_idx_t
#include "povu/refs/refs.hpp"	 // for pr, Ref, Refs

namespace povu::bidirected
{
inline constexpr std::string_view MODULE = "povu::bidirected";

using namespace povu::types::graph;
namespace pgt = povu::types::graph;

namespace lq = liteseq;

/**
 * undirected edge
 * stores the index of the vertex in the graph not the id
 */
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
	[[nodiscard]] pt::id_t id() const;
	[[nodiscard]] pt::u32 get_length() const; // length of the label
	[[nodiscard]] const std::string &get_label() const;
	// reverse complement of the label
	[[nodiscard]] std::string get_rc_label() const;
	[[nodiscard]] const std::set<pt::idx_t> &get_edges_l() const;
	[[nodiscard]] const std::set<pt::idx_t> &get_edges_r() const;

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

	/*
	  --------------------------
	  ref or hap related getters
	  --------------------------
	*/
	std::string get_sample_name(pt::id_t ref_id) const;
	std::string get_tag(pt::id_t ref_id) const;
	const pr::Ref &get_ref_by_id(pt::id_t ref_id) const;
	pr::Ref &get_ref_by_id_mut(pt::id_t ref_id);
	std::optional<pt::id_t> get_ref_id(std::string_view ref_tag) const;

	// bool has_refs() const;
	/** @brief if PanSN the prefix will be in sample */
	std::set<pt::id_t> get_shared_samples(pt::id_t ref_id) const;
	// sometimes the sample name is referred to as a prefix
	/** @brief the set of ref ids in the sample
	 *
	 * sometimes the sample name is referred to as a prefix
	 */
	std::set<pt::id_t>
	get_refs_in_sample(std::string_view sample_name) const;
	const lq::ref *get_ref_vec(pt::id_t ref_id) const;

	/** for GFA 1.1 returns the no of P lines in the graph */
	pt::u32 get_hap_count() const;

	const std::vector<pt::idx_t> &
	get_vertex_ref_idxs(pt::idx_t v_idx, pt::id_t ref_id) const;

	const std::vector<std::vector<pt::idx_t>> &
	get_vertex_refs(pt::idx_t v_id) const;

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
	void print_gfa(std::ostream &os) const;
};

using VG = VariationGraph;
// typedef VariationGraph VG;
} // namespace povu::bidirected

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace bd = povu::bidirected;

#endif
