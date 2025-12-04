#ifndef POVU_GENOMICS_ALLELE_HPP
#define POVU_GENOMICS_ALLELE_HPP

// #include <cstdlib> // for exit, EXIT_FAILURE
// #include <functional>
#include <liteseq/refs.h>	  // for ref_walk
#include <liteseq/types.h>	  // for strand
#include <map>			  // for map
#include <set>			  // for set, operator!=
#include <string>		  // for basic_string, string
#include <string_view>		  // for string_view
#include <utility>		  // for move, pair
#include <vector>		  // for vector
				  //
#include "povu/common/compat.hpp" // for contains, pv_cmp
#include "povu/common/constants.hpp"
#include "povu/common/core.hpp" // for pt, idx_t, id_t, op_t
// #include "povu/common/log.hpp"	  // for ERR
#include "povu/common/utils.hpp"
#include "povu/graph/bidirected.hpp" // for bd, VG
#include "povu/graph/pvst.hpp"	     // for VertexBase
#include "povu/graph/types.hpp"	     // for or_e, id_or_t, walk_t
// #include "povu/overlay/overlay.hpp"
#include "povu/variation/rov.hpp" // for RoV

namespace povu::genomics::allele
{
inline constexpr std::string_view MODULE = "povu::genomics::allele";

namespace lq = liteseq;
namespace pgt = povu::types::graph;
namespace pvst = povu::pvst;

// idx_in_hap, v_id, or_e
using extended_step = std::tuple<pt::u32, pt::u32, pgt::or_e>;
using extended_walk = std::vector<extended_step>;

using broad_step = extended_step;
using broad_walk = extended_walk; // a lap?
using ext_lap = std::vector<extended_step>;
using race = std::vector<broad_walk>;

struct hap_slice {
	const lq::ref_walk *ref_w;
	pt::idx_t ref_idx;
	pt::idx_t ref_start_idx;
	pt::idx_t len;

	// --------------
	// constructor(s)
	// --------------
	static hap_slice create_null()
	{
		return hap_slice{nullptr, pc::INVALID_IDX, pc::INVALID_IDX,
				 pc::INVALID_IDX};
	}

	hap_slice(const lq::ref_walk *ref_w_, pt::idx_t ref_idx_,
		  pt::idx_t ref_start_idx_, pt::idx_t len_)
	    : ref_w(ref_w_), ref_idx(ref_idx_), ref_start_idx(ref_start_idx_),
	      len(len_)
	{}

	hap_slice() = delete;

	// ---------
	// getter(s)
	// ---------

	[[nodiscard]]
	bool is_null() const
	{
		return this->ref_w == nullptr ||
		       this->ref_idx == pc::INVALID_IDX ||
		       this->ref_start_idx == pc::INVALID_IDX ||
		       this->len == pc::INVALID_IDX;
	}

	[[nodiscard]]
	bd::id_or_t get_step(pt::idx_t idx) const
	{
		pt::idx_t ref_v_id = ref_w->v_ids[idx];
		pgt::or_e ref_o = ref_w->strands[idx] == lq::strand::STRAND_FWD
					  ? pgt::or_e::forward
					  : pgt::or_e::reverse;

		return {ref_v_id, ref_o};
	}

	[[nodiscard]]
	std::string dbg_str() const
	{
		std::string at_str = "";

		pt::u32 i = ref_start_idx;
		pt::u32 N = ref_start_idx + len;

		for (; i < N; i++)
			at_str += this->get_step(i).as_str();

		return at_str;
	}

	[[nodiscard]]
	pgt::walk_t to_walk() const
	{
		pgt::walk_t w;
		pt::u32 i = ref_start_idx;
		pt::u32 N = ref_start_idx + len;

		for (; i < N; i++)
			w.emplace_back(this->get_step(i));

		return w;
	}

	[[nodiscard]]
	std::string get_step_label(const pgt::step_t &s, const bd::VG &g) const
	{
		auto [v_id, o] = s;
		return (o == pgt::or_e::forward)
			       ? g.get_vertex_by_id(v_id).get_label()
			       : g.get_vertex_by_id(v_id).get_rc_label();
	}

	// void append_step_label(const pgt::step_t &s, const bd::VG &g,
	//		       std::string &dna_str) const
	// {
	//	auto [v_id, o] = s;
	//	dna_str += (o == pgt::or_e::forward)
	//			   ? g.get_vertex_by_id(v_id).get_label()
	//			   : g.get_vertex_by_id(v_id).get_rc_label();
	// }

	[[nodiscard]]
	std::string as_str(pvr::var_type_e vt) const
	{
		std::string at_str = "";

		pt::u32 i = ref_start_idx;
		pt::u32 N = ref_start_idx + len;

		switch (vt) {
		case pvr::var_type_e::sub:
			i++;
			N--;
			break;
		case pvr::var_type_e::ins:
		case pvr::var_type_e::del:
			N--;
			break;
		case pvr::var_type_e::inv:
			break;
		}

		for (; i < N; i++)
			at_str += this->get_step(i).as_str();

		return at_str;
	}

	[[nodiscard]]
	pt::u32 comp_pos(pvr::var_type_e vt) const
	{
		pt::idx_t locus = ref_w->loci[ref_start_idx + 1];

		switch (vt) {
		case pvr::var_type_e::del:
		case pvr::var_type_e::ins:
			return locus - 1;
		case pvr::var_type_e::sub:
		case pvr::var_type_e::inv:
			return locus;
		}

		ERR("Unknown variant type");

		return pc::INVALID_IDX;
	}

	[[nodiscard]]
	std::string as_dna_str(const bd::VG &g, pvr::var_type_e vt) const
	{
		std::string dna_str = "";

		pt::u32 i = ref_start_idx;
		pt::u32 N = ref_start_idx + len;

		switch (vt) {
		case pvr::var_type_e::inv:
		case pvr::var_type_e::sub:
			break;
		case pvr::var_type_e::ins:
		case pvr::var_type_e::del:;
			const pgt::step_t &s = this->get_step(i);
			dna_str += get_step_label(s, g).back();
			break;
		}

		if (vt != pvr::var_type_e::inv) {
			i++;
			N--;
		}

		for (; i < N; i++)
			dna_str += get_step_label(this->get_step(i), g);

		return dna_str;
	}
};

// used for untangle
// TODO: rename to allele traversal
struct at_itn {
	std::vector<pgt::walk_t> it_;

	// --------------
	// constructor(s)
	// --------------
	at_itn() : it_()
	{}

	at_itn(std::vector<pgt::walk_t> &&itn) : it_(itn)
	{}

	// ---------
	// getter(s)
	// ---------
	[[nodiscard]]
	pt::idx_t at_count() const
	{
		return this->it_.size();
	}

	[[nodiscard]]
	const std::vector<pgt::walk_t> &get_ats() const
	{
		return this->it_;
	}

	[[nodiscard]]
	const pgt::walk_t &get_at(pt::idx_t at_idx) const
	{
		return this->it_[at_idx];
	}

	std::string to_string()
	{
		std::string res = "";
		for (const pgt::walk_t &w : this->it_)
			res += pv_cmp::format("{{ {} }}", pgt::to_string(w));

		return res;
	}
};

struct rov_boundaries {
private:
	pgt::id_or_t l_;
	pgt::id_or_t r_;

public:
	rov_boundaries() = delete;

	rov_boundaries(pgt::id_or_t l, pgt::id_or_t r) : l_(l), r_(r)
	{}

	[[nodiscard]]
	pt::op_t<pgt::id_or_t> get_bounds() const
	{
		return pt::op_t<pgt::id_or_t>{this->l_, this->r_};
	}

	[[nodiscard]]
	std::string to_string() const
	{
		return pv_cmp::format("{}{}", this->l_.as_str(),
				      this->r_.as_str());
	}

	static rov_boundaries create_null()
	{
		return rov_boundaries{
			pgt::id_or_t{pc::INVALID_IDX, pgt::or_e::forward},
			pgt::id_or_t{pc::INVALID_IDX, pgt::or_e::forward}};
	}

	bool friend operator<(const rov_boundaries &lhs,
			      const rov_boundaries &rhs)
	{
		return std::make_pair(lhs.l_, lhs.r_) <
		       std::make_pair(rhs.l_, rhs.r_);
	}

	std::ostream friend &operator<<(std::ostream &os,
					const rov_boundaries &cxt)
	{
		return os << cxt.to_string();
	}
};

using walk_to_alts_map = std::map<pgt::walk_t, std::vector<hap_slice>>;

struct alt_set {
private:
	walk_to_alts_map ins;
	walk_to_alts_map dels;
	walk_to_alts_map subs;

public:
	alt_set() = default;

	void add_del(hap_slice &&s)
	{
		this->dels[s.to_walk()].emplace_back(s);
	}

	[[nodiscard]]
	const walk_to_alts_map &get_ins() const
	{
		return this->ins;
	}

	[[nodiscard]]
	const walk_to_alts_map &get_dels() const
	{
		return this->dels;
	}

	[[nodiscard]]
	const walk_to_alts_map &get_subs() const
	{
		return this->subs;
	}

	void add_ins(hap_slice &&s)
	{
		this->ins[s.to_walk()].emplace_back(s);
	}

	void add_sub(hap_slice &&s)
	{
		this->subs[s.to_walk()].emplace_back(s);
	}

	void print_alt(std::ostream &os, pvr::var_type_e vt) const
	{
		const walk_to_alts_map *alt = nullptr;

		switch (vt) {
		case pvr::var_type_e::sub:
			alt = &this->subs;
			break;
		case pvr::var_type_e::ins:
			alt = &this->ins;
			break;
		case pvr::var_type_e::del:
			alt = &this->dels;
			break;
		case pvr::var_type_e::inv:
			return;
		}

		for (const auto &[w, alts] : *alt) {
			os << "Walk: " << pgt::to_string(w) << "\t Alts Haps: ";
			for (auto it = alts.begin(); it != alts.end(); ++it) {
				os << it->ref_idx;
				if (std::next(it) != alts.end())
					os << ", ";
			}
			os << "\n";
		}
	}
};

struct minimal_rov {
private:
	rov_boundaries cxt_ = rov_boundaries::create_null(); // context
	hap_slice ref_as_ = hap_slice::create_null(); // reference allele slice
	std::set<pt::u32> haps_matching_ref;
	std::set<pt::u32> alt_haps;

	alt_set alts_; // alternative allele slices

public:
	// --------------
	// constructor(s)
	// --------------

	minimal_rov() = delete;

	minimal_rov(const rov_boundaries cxt, hap_slice &&ref_as)
	    : cxt_(cxt), ref_as_(ref_as), alts_()
	{}

	// ---------
	// getter(s)
	// ---------

	[[nodiscard]] bool ref_as_is_null() const
	{
		return this->ref_as_.is_null();
	}

	[[nodiscard]] hap_slice get_ref_as() const
	{
		return this->ref_as_;
	}

	[[nodiscard]] rov_boundaries get_context() const
	{
		return this->cxt_;
	}

	[[nodiscard]] pt::u32 get_ref_len() const
	{
		return this->get_ref_as().len;
	}

	[[nodiscard]]
	const std::set<pt::u32> &get_haps_matching_ref() const
	{
		return this->haps_matching_ref;
	}

	[[nodiscard]]
	const std::set<pt::u32> &get_alt_haps() const
	{
		return this->alt_haps;
	}

	[[nodiscard]]
	pt::u32 get_ref_at_ref_count() const
	{
		return this->haps_matching_ref.size();
	}

	void add_haps_match_ref(pt::u32 h_idx)
	{
		this->haps_matching_ref.insert(h_idx);
	}

	[[nodiscard]]
	const walk_to_alts_map &get_subs() const
	{
		return this->alts_.get_subs();
	}

	[[nodiscard]]
	const walk_to_alts_map &get_ins() const
	{
		return this->alts_.get_ins();
	}

	[[nodiscard]]
	const walk_to_alts_map &get_dels() const
	{
		return this->alts_.get_dels();
	}

	// rename to add_alt_slice
	void add_alt(hap_slice &&alt_as)
	{
		this->alt_haps.insert(alt_as.ref_idx);

		if (alt_as.len == this->get_ref_len()) {
			this->alts_.add_sub(std::move(alt_as));
		}
		else if (alt_as.len > this->get_ref_len()) { // insertion
			this->alts_.add_ins(std::move(alt_as));
		}
		else { // deletion
			this->alts_.add_del(std::move(alt_as));
		}
	}

	void print(std::ostream &os) const
	{
		os << "minimal Rov:\n";
		os << this->get_context() << "\n";
		os << "Ref " << this->get_ref_as().dbg_str();
		os << "\t alts: ";
		this->alts_.print_alt(os, pvr::var_type_e::sub);
		this->alts_.print_alt(os, pvr::var_type_e::ins);
		this->alts_.print_alt(os, pvr::var_type_e::del);
	}
};

using cxt_to_min_rov_map = std::map<rov_boundaries, minimal_rov>;

struct trek {
private:
	// pointer to the RoV from which the expedition is made
	const pvr::RoV *rov_;

	// hap idx to context to minimal Rov map
	// std::vector<cxt_to_min_rov_map *> data;
	std::map<pt::u32, cxt_to_min_rov_map> data_;

	// key is the ref and value is the set of haps that do not cover the RoV
	std::map<pt::u32, std::set<pt::u32>> no_cov; // do not cover at all

	// key is the ref and value is the set of haps that match the ref
	std::map<pt::u32, std::set<pt::u32>> matches_ref;

	bool tangled_{false}; // is true when tangling exists
	const pt::u32 HAP_COUNT{pc::INVALID_IDX};

	// ----------------------
	// private constructor(s)
	// ----------------------
	// trek(const pvr::RoV *rov, std::vector<cxt_to_min_rov_map *> d,
	//      pt::u32 hap_count, bool tangled = false)
	//     : rov_(rov), data(std::move(d)), tangled_(tangled),
	//       HAP_COUNT(hap_count)
	// {}

	trek(const pvr::RoV *rov, pt::u32 hap_count, bool tangled = false)
	    : rov_(rov), tangled_(tangled), HAP_COUNT(hap_count)
	{}

public:
	// --------------
	// constructor(s)
	// --------------
	trek() = delete;

	[[nodiscard]]
	static trek create_new(const pvr::RoV *rov, pt::u32 hap_count,
			       bool tangled)
	{
		return {rov, hap_count, tangled};
	}

	// ---------
	// getter(s)
	// ---------

	[[nodiscard]]
	const pvst::VertexBase *get_pvst_vtx_const_ptr() const
	{
		return this->rov_->get_pvst_vtx();
	}

	[[nodiscard]]
	cxt_to_min_rov_map &get_min_rov(pt::u32 h_idx)
	{
		return this->data_[h_idx];
	}

	// [[nodiscard]]
	// bool has_context(pt::u32 ref_h_idx, const rov_boundaries cxt)
	// {
	//	cxt_to_min_rov_map *d = data[ref_h_idx];
	//	return pv_cmp::contains(*d, cxt);
	// };

	[[nodiscard]]
	bool has_data() const
	{
		return !this->data_.empty();
	}

	// [[nodiscard]]
	// bool has_data() const
	// {
	//	for (auto d : data)
	//		if (d != nullptr)
	//			return true;

	//	return false;
	// }

	[[nodiscard]]
	bool is_tangled() const
	{
		return this->tangled_;
	}

	[[nodiscard]]
	pt::u32 get_hap_count() const
	{
		return this->HAP_COUNT;
	}

	// TODO: replace with a set member of the struct
	[[nodiscard]]
	std::set<pt::u32> get_ref_haps() const
	{
		std::set<pt::u32> ref_haps;
		for (const auto &[ref_h_idx, _] : this->data_)
			ref_haps.insert(ref_h_idx);

		// for (pt::u32 ref_h_idx{}; ref_h_idx < HAP_COUNT; ref_h_idx++)
		//	if (data[ref_h_idx] != nullptr)
		//		ref_haps.insert(ref_h_idx);

		return ref_haps;
	}

	[[nodiscard]]
	const cxt_to_min_rov_map &get_ref_recs(pt::u32 ref_h_idx) const
	{
		return this->data_.at(ref_h_idx);
	}

	[[nodiscard]]
	cxt_to_min_rov_map &get_ref_recs_mut(pt::u32 ref_h_idx)
	{
		return this->data_.at(ref_h_idx);
	}

	[[nodiscard]]
	const std::set<pt::u32> get_match_ref(pt::u32 ref_h_idx) const
	{
		// TODO: handle key not existing
		return this->matches_ref.at(ref_h_idx);
	}

	[[nodiscard]]
	const std::set<pt::u32> get_no_cov(pt::u32 ref_h_idx) const
	{
		// TODO: handle key not existing
		return this->no_cov.at(ref_h_idx);
	}

	// ---------
	// setter(s)
	// ---------

	// void init_ref_idx(pt::u32 ref_h_idx)
	// {
	//	// TODO: [A] memory leak
	//	if (data[ref_h_idx] == nullptr)
	//		data[ref_h_idx] = new cxt_to_min_rov_map();
	// }

	void add_no_cov(pt::u32 ref_h_idx, pt::u32 h_idx)
	{
		this->no_cov[ref_h_idx].insert(h_idx);
	}

	void add_match_ref(pt::u32 ref_h_idx, pt::u32 h_idx)
	{
		this->matches_ref[ref_h_idx].insert(h_idx);
	}

	void set_tangled(bool t)
	{
		this->tangled_ = t;
	}

	// ---------
	// debug
	// ---------

	void dbg_print(std::ostream &os) const
	{
		for (const auto &[ref_h_idx, d] : this->data_) {
			os << "Ref hap idx " << ref_h_idx << "\n";

			os << "No coverage: ";
			if (pv_cmp::contains(no_cov, ref_h_idx)) {
				os << "{";
				std::cerr << pu::concat_with(
					no_cov.at(ref_h_idx), ',');
				std::cerr << "}\n";
			}
			else {
				os << "none\n";
			}

			os << "Records:\n";
			for (const auto &[cxt, min_rov] : d) {
				min_rov.print(os);
				os << "\n";
			}
		}

		// for (pt::u32 ref_h_idx{}; ref_h_idx < data.size();
		//      ref_h_idx++) {
		//	cxt_to_min_rov_map *d = data[ref_h_idx];
		//	if (d == nullptr)
		//		continue;

		//	os << "Ref hap idx " << ref_h_idx << "\n";

		//	os << "No coverage: ";
		//	if (pv_cmp::contains(no_cov, ref_h_idx)) {
		//		os << "{";
		//		std::cerr << pu::concat_with(
		//			no_cov.at(ref_h_idx), ',');
		//		std::cerr << "}\n";
		//	}
		//	else {
		//		os << "none\n";
		//	}

		//	os << "Records:\n";
		//	for (const auto &[cxt, min_rov] : *d) {
		//		min_rov.print(os);
		//		os << "\n";
		//	}
		// }
	}
};

// precense absence matrix
struct depth_matrix {
private:
	std::vector<pt::u32> data;

	pt::u32 I; // haps (rows)
	pt::u32 J; // vertices in BFS tree (cols)

	pt::u32 max_depth{0};

	std::map<pt::u32, pt::u32> hap_idx_to_loop_no;

	std::vector<pt::u32> sorted_vertices;

	bool is_tangled{false};

public:
	// --------------
	// constructor(s)
	// --------------
	// store in row major order

	depth_matrix() = delete;

	[[nodiscard]]
	depth_matrix(pt::u32 i, pt::u32 j)
	    : data(i * j, 0), I(i), J(j)
	{}

	// use this for debugging especially
	[[nodiscard]]
	depth_matrix(pt::u32 i, pt::u32 j, const std::vector<pt::u32> &sv)
	    : data(i * j, 0), I(i), J(j), sorted_vertices(sv)
	{}

	// ---------
	// getter(s)
	// ---------

	[[nodiscard]]
	const std::vector<pt::u32> &get_header() const
	{
		return this->sorted_vertices;
	}

	[[nodiscard]]
	std::vector<pt::u32> get_col_data(pt::u32 j) const
	{
		std::vector<pt::u32> col_data(I, 0);
		for (pt::u32 i{}; i < I; i++)
			col_data[i] = data[i * J + j];

		return col_data;
	}

	[[nodiscard]]
	std::vector<pt::u32> get_row_data(pt::u32 i) const
	{
		std::vector<pt::u32> row_data;
		row_data.reserve(J);
		for (pt::u32 j{}; j < J; j++)
			row_data.push_back(data[i * J + j]);

		return row_data;
	}

	void set_loop_no(pt::u32 h_idx, pt::u32 loop_no)
	{
		this->hap_idx_to_loop_no[h_idx] = loop_no;
	}

	[[nodiscard]]
	pt::u32 get_loop_no(pt::u32 h_idx) const
	{
		if (!this->is_tangled)
			return 0;

		if (pv_cmp::contains(this->hap_idx_to_loop_no, h_idx))
			return this->hap_idx_to_loop_no.at(h_idx);

		return pc::INVALID_IDX;
	}

	[[nodiscard]]
	pt::u32 row_count() const
	{
		return I;
	}

	[[nodiscard]]
	pt::u32 col_count() const
	{
		return J;
	}

	bool tangled()
	{
		return is_tangled;
	}

	[[nodiscard]]
	pt::u32 get_depth(pt::u32 i, pt::u32 j) const
	{
		return this->get_data()[i * J + j];
	}

	[[nodiscard]]
	pt::u32 get_max_depth() const
	{
		return this->max_depth;
	}

	// -----------
	// modifier(s)
	// -----------

	void set_tangled(bool t)
	{
		this->is_tangled = t;
	}

	void fill(const bd::VG &g, const pvr::RoV &rov)
	{
		// fill column-wise
		// row wise is more cache friendly, but needs
		// a method in bd::VG to get haplotypes per haplotype
		for (pt::u32 h_idx{}; h_idx < I; h_idx++) {
			for (pt::u32 j{}; j < J; j++) {
				pt::u32 v_id = rov.get_sorted_vertex(j);

				pt::u32 v_idx = g.v_id_to_idx(v_id);
				const std::vector<pt::idx_t> &ref_idxs =
					g.get_vertex_ref_idxs(v_idx, h_idx);
				pt::u32 depth = ref_idxs.size();

				if (depth == 0)
					continue;

				if (depth == 1) {
					const lq::ref_walk *h_w =
						g.get_ref_vec(h_idx)
							->walk;	 // the hap walk
					pt::u32 k = ref_idxs[0]; // index in the
								 // hap walk
					pgt::or_e orn =
						h_w->strands[k] ==
								lq::strand::
									STRAND_FWD
							? pgt::or_e::forward
							: pgt::or_e::reverse;

					switch (orn) {
					case pgt::or_e::forward:
						this->set_data(
							h_idx,
							rov.get_sorted_pos(
								v_id),
							1);
						break;
					case pgt::or_e::reverse:
						this->set_data(
							h_idx,
							rov.get_sorted_pos(
								v_id),
							2);
						break;
					}
					continue;
				}

				if (depth > 1 && (j == 0 || j == J - 1))
					this->is_tangled = true;

				if (depth > this->max_depth)
					this->max_depth = depth;

				this->set_data(h_idx, j, depth);
			}
		}
	}

	[[nodiscard]]
	const std::vector<pt::u32> get_data() const
	{
		return this->data;
	}

	[[nodiscard]]
	std::vector<pt::u32> &get_data_mut()
	{
		return this->data;
	}

	void set_data(pt::u32 i, pt::u32 j, pt::u32 d)
	{
		this->data[i * J + j] = d;
	}

	// other
	void print(std::ostream &os) const
	{
		os << "is tangled " << (this->is_tangled ? "true" : "false")
		   << "\n";

		// print hap_idx_to_loop_no map
		os << "#HapIdxToLoopNo\n";
		for (const auto &[h_idx, loop_no] : this->hap_idx_to_loop_no)
			os << pv_cmp::format("[hap {}, loop {}] ", h_idx,
					     loop_no);

		os << "\n";

		// print header
		os << "-\t" << pu::concat_with(this->sorted_vertices, '\t');
		os << "\n";

		// print rows
		for (pt::u32 i{}; i < I; i++) {
			os << i << "\t"
			   << pu::concat_with(get_row_data(i), '\t');
			os << "\n";
		}
	}
};

} // namespace povu::genomics::allele

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pga = povu::genomics::allele;

#endif // POVU_GENOMICS_ALLELE_HPP
