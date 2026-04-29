#ifndef BIDIRECTED_HPP
#define BIDIRECTED_HPP

#include <cstddef>     // for size_t
#include <iostream>    // for ostream
#include <optional>    // for optional
#include <set>	       // for set
#include <string>      // for string, basic_string
#include <string_view> // for string_view
#include <vector>      // for vector

#include <liteseq/gfa.h>	 // for gfa_free, gfa_props
#include <liteseq/refs.h>	 // for ref
#include <meza/owned/matrix.hpp> // for dense_matrix2d

#include <quilt/types.hpp> // for qt

#include "povu/common/utils.hpp" // for pu, TwoWayMap
#include "povu/graph/types.hpp"	 // for v_end_e, side_n_id_t, side_n_idx_t
#include "povu/refs/refs.hpp"	 // for pr, Ref, Refs

namespace oza::bidirected
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
	qt::idx_t v1_idx_;
	pgt::v_end_e v1_end_;
	qt::idx_t v2_idx_;
	pgt::v_end_e v2_end_;

public:
	// --------------
	// constructor(s)
	// --------------
	Edge(qt::idx_t v1_idx, pgt::v_end_e v1_end, qt::idx_t v2_idx,
	     pgt::v_end_e v2_end);

	// ---------
	// getter(s)
	// ---------
	qt::idx_t get_v1_idx() const;
	qt::idx_t &get_v1_idx_mut();
	pgt::v_end_e get_v1_end() const;
	qt::idx_t get_v2_idx() const;
	qt::idx_t &get_v2_idx_mut();
	pgt::v_end_e get_v2_end() const;
	pgt::side_n_idx_t get_other_vtx(qt::idx_t v_id,
					pgt::v_end_e v_end) const;
	pgt::side_n_idx_t
	get_other_vtx(qt::idx_t v_id) const; // if you don't care for self loops
};

class Vertex
{
	qt::id_t v_id_;	    // this is the sequence name in the GFA file. Maybe
			    // Should support strings.
	std::string label_; // or sequence

	// indexes to the edge vector in Graph
	std::set<qt::idx_t> e_l;
	std::set<qt::idx_t> e_r;

public:
	// --------------
	// constructor(s)
	// --------------
	Vertex(qt::id_t v_id, const std::string &label = "");

	// ---------
	// getter(s)
	// ---------
	[[nodiscard]] qt::id_t id() const;
	[[nodiscard]] qt::u32 get_length() const; // length of the label
	[[nodiscard]] const std::string &get_label() const;
	// reverse complement of the label
	[[nodiscard]] std::string get_rc_label() const;
	[[nodiscard]] const std::set<qt::idx_t> &get_edges_l() const;
	[[nodiscard]] const std::set<qt::idx_t> &get_edges_r() const;

	// ---------
	// setter(s)
	// ---------
	void add_edge_l(qt::idx_t e_idx);
	void add_edge_r(qt::idx_t e_idx);
};

class VariationGraph
{
	std::vector<Vertex> vertices;
	std::set<pgt::side_n_id_t> tips_; // the set of side and id of the tips
	pu::TwoWayMap<std::size_t, std::size_t>
		v_id_to_idx_; // TODO: reserve size

	std::vector<Edge> edges;

	// vtx2sm -> vertex to step matrix
	// i is the vertex index, j is the ref index, k is the step index
	meza::matrix::dense_matrix2d<std::vector<qt::idx_t>> vtx2sm;

	pr::Refs refs_;

	lq::gfa_props *gfa;

public:
	// --------------
	// constructor(s)
	// --------------
	VariationGraph(qt::idx_t vtx_count, qt::idx_t edge_count,
		       qt::idx_t ref_count);
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
	qt::idx_t v_id_to_idx(qt::id_t v_id) const;
	qt::id_t v_idx_to_id(qt::idx_t v_idx) const;

	qt::idx_t vtx_count() const;
	qt::idx_t edge_count() const;
	const std::set<pgt::side_n_id_t> &tips() const;
	const Edge &get_edge(qt::idx_t e_idx) const;
	Edge &get_edge_mut(qt::idx_t e_idx);
	// TODO replace vertex with v?
	const Vertex &get_vertex_by_idx(qt::idx_t v_idx) const;
	const Vertex &get_vertex_by_id(qt::id_t v_id) const;
	Vertex &get_vertex_mut_by_id(qt::id_t v_id);

	/*
	  --------------------------
	  ref or hap related getters
	  --------------------------
	*/
	std::string get_sample_name(qt::id_t ref_id) const;
	std::string get_tag(qt::id_t ref_id) const;
	const pr::Ref &get_ref_by_id(qt::id_t ref_id) const;
	pr::Ref &get_ref_by_id_mut(qt::id_t ref_id);
	std::optional<qt::id_t> get_ref_id(std::string_view ref_tag) const;

	// bool has_refs() const;
	/** @brief if PanSN the prefix will be in sample */
	std::set<qt::id_t> get_shared_samples(qt::id_t ref_id) const;
	// sometimes the sample name is referred to as a prefix
	/** @brief the set of ref ids in the sample
	 *
	 * sometimes the sample name is referred to as a prefix
	 */
	std::set<qt::id_t>
	get_refs_in_sample(std::string_view sample_name) const;
	const lq::ref *get_ref_vec(qt::id_t ref_id) const;

	/** for GFA 1.1 returns the no of P lines in the graph */
	qt::u32 get_hap_count() const;

	const std::vector<qt::idx_t> &
	get_vertex_ref_idxs(qt::idx_t v_idx, qt::id_t ref_id) const;

	qt::idx_t get_ploidy(const std::string &sample_name) const;
	qt::idx_t get_ploidy_id(const std::string &sample_name,
				qt::u32 ploidy_idx) const;

	const std::vector<std::string> &get_genotype_col_names() const;
	std::vector<std::vector<std::string>> get_blank_genotype_cols() const;
	const pr::gt_col_meta &get_gt_col_meta(qt::id_t ref_id) const;

	// ---------
	// setter(s)
	// ---------
	void add_tip(qt::id_t v_id, pgt::v_end_e end);
	// returns the index (v_idx) of the added vertex
	qt::idx_t add_vertex(qt::id_t v_id, const std::string &label);
	// returns the index (e_idx) of the added edge
	qt::idx_t add_edge(qt::id_t v1_id, pgt::v_end_e v1_end, qt::id_t v2_id,
			   pgt::v_end_e v2_end);
	void set_refs_meta(lq::ref **refs, qt::idx_t ref_count);
	// qt::id_t add_ref(const std::string &label, char delim);
	void shrink_to_fit();
	void set_vtx_ref_idx(qt::id_t v_id, qt::id_t ref_id,
			     qt::idx_t step_idx);

	// -----
	// other
	// -----
	void summary(bool print_tips) const;
	void print_dot(std::ostream &os) const;
	void print_gfa(std::ostream &os) const;
};

using VG = VariationGraph;
} // namespace oza::bidirected

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace bd = oza::bidirected;

#endif
