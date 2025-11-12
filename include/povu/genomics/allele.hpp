#ifndef POVU_GENOMICS_ALLELE_HPP
#define POVU_GENOMICS_ALLELE_HPP

#include <cstdlib>		     // for exit, EXIT_FAILURE
#include <liteseq/refs.h>	     // for ref_walk
#include <liteseq/types.h>	     // for strand
#include <map>			     // for map
#include <set>			     // for set, operator!=
#include <string>		     // for basic_string, string
#include <string_view>		     // for string_view
#include <utility>		     // for move, pair
#include <vector>		     // for vector
				     //
#include "povu/common/compat.hpp"    // for contains, pv_cmp
#include "povu/common/core.hpp"	     // for pt, idx_t, id_t, op_t
#include "povu/common/log.hpp"	     // for ERR
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

struct ref_needle {
	pt::u32 r_idx;
	pt::u32 start;
	pt::u32 limit_left{pc::INVALID_IDX};
	pt::u32 limit_right{pc::INVALID_IDX};

	ref_needle(pt::u32 r_idx_, pt::u32 start_)
	    : r_idx(r_idx_), start(start_)
	{}

	void set_limits(pt::u32 left, pt::u32 right)
	{
		limit_left = left;
		limit_right = right;
	}
};

struct sub_inv {
	std::vector<ref_needle> rev_needles;
	std::vector<ref_needle> fwd_needles;
	pt::u32 len;
	const pvr::RoV *rov_;
};

struct allele_slice_t {
	const pgt::walk_t *walk;
	pt::idx_t walk_idx;
	pt::idx_t walk_start_idx;

	const lq::ref_walk *ref_w;
	pt::idx_t ref_idx;
	pt::idx_t ref_start_idx;

	pt::idx_t len; // total step count in the itinerary
	ptg::or_e slice_or;
	pvr::var_type_e vt;

	// ---------
	// getter(s)
	// ---------
	[[nodiscard]]
	pt::idx_t step_count() const
	{
		return this->len;
	}

	[[nodiscard]]
	ptg::or_e get_or() const
	{
		return this->slice_or;
	}

	[[nodiscard]]
	bd::id_or_t get_walk_step(pt::idx_t idx) const
	{
		return this->walk->at(idx);
	}

	[[nodiscard]]
	pt::idx_t get_locus(pt::idx_t idx) const
	{
		return this->ref_w->loci[idx];
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
	std::string as_str(pvr::var_type_e variant_type) const
	{
		// TODO: pass variant type as param this method is prone to bugs
		bool is_fwd = this->get_or() == pgt::or_e::forward;
		std::string at_str = "";

		pt::u32 i = is_fwd ? ref_start_idx : ref_start_idx - len + 1;
		pt::u32 N = is_fwd ? ref_start_idx + len : ref_start_idx + 1;

		switch (variant_type) {
		case pvr::var_type_e::sub:
			i++;
			N--;
			break;
		case pvr::var_type_e::ins:
		case pvr::var_type_e::del:
			N--;
			break;
		case pvr::var_type_e::und: // undefined
			ERR("Undefined variant type in allele_slice_t::as_str");
			std::exit(EXIT_FAILURE);
		}

		for (; i < N; i++)
			at_str += this->get_step(i).as_str();

		return at_str;
	}
};

bool ref_eq(const allele_slice_t &lhs, const allele_slice_t &rhs);
bool operator==(const allele_slice_t &lhs, const allele_slice_t &rhs);
bool operator!=(const allele_slice_t &lhs, const allele_slice_t &rhs);

/**
 * Ref Itinerary or just Itinerary
 * an interrupted
 * sequence of looped walks in a RoV for a given ref
 * useful for repeats
 */
struct itn_t {
	std::vector<allele_slice_t> it_;

	// --------------
	// constructor(s)
	// --------------
	itn_t() : it_()
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
	const std::vector<allele_slice_t> &get_ats() const
	{
		return this->it_;
	}

	[[nodiscard]]
	const allele_slice_t &get_at(pt::idx_t at_idx) const
	{
		return this->it_[at_idx];
	}

	// ---------
	// setter(s)
	// ---------
	void append_at(allele_slice_t &&s)
	{
		this->it_.emplace_back(s);
	}

	void append_at_sorted(allele_slice_t &&s)
	{
		if (this->it_.empty()) {
			this->it_.emplace_back(s);
			return;
		}

		// sort in ascending order based on ref_start_idx
		pt::idx_t insert_idx = 0;
		for (const auto &at : this->it_) {
			if (s.ref_start_idx < at.ref_start_idx) {
				break;
			}
			insert_idx++;
		}
		this->it_.insert(this->it_.begin() + insert_idx, s);
	}

	void sort()
	{
		std::sort(this->it_.begin(), this->it_.end(),
			  [](const allele_slice_t &a, const allele_slice_t &b)
			  { return a.ref_start_idx < b.ref_start_idx; });
	}
};

/**
 * Exp, short for expedition: a journey undertaken by a group of people with
 * a particular purpose, especially that of exploration, research, or war.
 *
 * A collection of itineraries for each reference in a region of variation
 * map of ref_id to the itn of the ref in a RoV
 */
class Exp
{
	// pointer to the RoV from which the expedition is made
	const pvr::RoV *rov_;

	// map of ref_id to the itinerary (set of walks) of the ref in a RoV
	// when tangled, a ref can have multiple walks in a RoV
	std::map<pt::id_t, itn_t> ref_itns_;

	// walk idx to ref idxs that take the walk
	std::map<pt::idx_t, std::set<pt::idx_t>> walk_idx_to_ref_idxs_;

	// alignment between two refs
	std::map<pt::op_t<pt::id_t>, std::string> aln_;

	// is true when tangling exists.
	// tangling exists when a walk traverses an RoV more than once
	bool is_tangled_{false};

public:
	// ---------------------
	// public constructor(s)
	// ---------------------

	Exp() : rov_(nullptr), ref_itns_(), walk_idx_to_ref_idxs_(), aln_()
	{}

	Exp(const pvr::RoV *rov)
	    : rov_(rov), ref_itns_(), walk_idx_to_ref_idxs_(), aln_()
	{
		if (this->rov_ == nullptr) {
			ERR("RoV pointer is null");
			std::exit(EXIT_FAILURE);
		}
	}

	// // disable copy constructor
	// Exp(const Exp &other) = delete;
	// Exp &operator=(const Exp &other) = delete;

	// // move constructor
	// Exp(Exp &&other) noexcept = default;
	// Exp &operator=(Exp &&other) noexcept = default;

	// ~Exp() = default;

	// ---------
	// getter(s)
	// ---------
	[[nodiscard]]
	pt::idx_t ref_count() const
	{
		return this->ref_itns_.size();
	}

	[[nodiscard]]
	const pvst::VertexBase *get_pvst_vtx_const_ptr() const
	{
		return this->rov_->get_pvst_vtx();
	}

	[[nodiscard]]
	std::string id() const
	{
		return this->rov_->as_str();
	}

	[[nodiscard]]
	std::set<pt::id_t> get_walk_idxs() const
	{
		std::set<pt::id_t> walk_idxs;
		for (const auto &p : this->walk_idx_to_ref_idxs_) {
			walk_idxs.insert(p.first);
		}
		return walk_idxs;
	}

	[[nodiscard]]
	std::set<pt::id_t> get_ref_ids() const
	{
		std::set<pt::id_t> ref_ids;
		for (const auto &p : this->ref_itns_) {
			ref_ids.insert(p.first);
		}
		return ref_ids;
	}

	[[nodiscard]]
	const itn_t &get_itn(pt::id_t ref_id) const
	{
		return this->ref_itns_.at(ref_id);
	}

	[[nodiscard]]
	itn_t &get_itn_mut(pt::id_t ref_id)
	{
		return this->ref_itns_.at(ref_id);
	}

	[[nodiscard]]
	const pvr::RoV *get_rov() const
	{
		return this->rov_;
	}

	[[nodiscard]]
	pt::idx_t walk_count() const
	{
		return this->rov_->walk_count();
	}

	std::map<pt::idx_t, std::set<pt::idx_t>> &get_walk_to_ref_idxs_mut()
	{
		return this->walk_idx_to_ref_idxs_;
	}

	[[nodiscard]]
	const std::set<pt::idx_t> &
	get_ref_idxs_for_walk(const pt::idx_t walk_idx) const
	{
		return this->walk_idx_to_ref_idxs_.at(walk_idx);
	}

	[[nodiscard]]
	const std::map<pt::id_t, itn_t> &get_ref_itns() const
	{
		return this->ref_itns_;
	}

	[[nodiscard]]
	std::map<pt::id_t, itn_t> &get_ref_itns_mut()
	{
		return this->ref_itns_;
	}

	[[nodiscard]]
	const std::string &get_aln(pt::id_t ref_id1, pt::id_t ref_id2) const
	{
		return this->aln_.at(pt::op_t<pt::id_t>{ref_id1, ref_id2});
	}

	[[nodiscard]]
	bool has_aln(pt::id_t ref_id1, pt::id_t ref_id2) const
	{
		return pv_cmp::contains(this->aln_,
					pt::op_t<pt::id_t>{ref_id1, ref_id2});
	}

	[[nodiscard]]
	const std::map<pt::op_t<pt::id_t>, std::string> &get_alns() const
	{
		return this->aln_;
	}

	[[nodiscard]]
	bool is_tangled() const
	{
		return this->is_tangled_;
	}

	[[nodiscard]]
	bool has_ref(pt::id_t ref_id) const
	{
		return pv_cmp::contains(this->ref_itns_, ref_id);
	}

	// ---------
	// setter(s)
	// ---------
	void add_aln(pt::id_t ref_id1, pt::id_t ref_id2, std::string &&aln)
	{
		this->aln_[pt::op_t<pt::id_t>{ref_id1, ref_id2}] = aln;
	}

	void set_tangled(bool is_tangled)
	{
		this->is_tangled_ = is_tangled;
	}
};

// std::pair<std::vector<Exp>, std::vector<sub_inv>>
// comp_itineraries3(const bd::VG &g, const pvr::RoV &rov,
//		  const std::set<pt::id_t> &to_call_ref_ids);

// std::pair<std::vector<Exp>, std::vector<pos::pin_cushion>>
// comp_itineraries3(const bd::VG &g, const pvr::RoV &rov,
//		  const std::set<pt::id_t> &to_call_ref_ids);

// std::vector<Exp> comp_itineraries2(const bd::VG &g, const pvr::RoV &rov);
// void comp_itineraries(const bd::VG &g, Exp &exp);

} // namespace povu::genomics::allele

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pga = povu::genomics::allele;

#endif // POVU_GENOMICS_ALLELE_HPP
