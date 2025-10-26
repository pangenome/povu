#ifndef POVU_GENOMICS_ALLELE_HPP
#define POVU_GENOMICS_ALLELE_HPP

#include <cstdlib>     // for exit, EXIT_FAILURE
#include <map>	       // for map
#include <set>	       // for set, operator!=
#include <string>      // for basic_string, string
#include <string_view> // for string_view
#include <utility>     // for move, pair
#include <vector>      // for vector

#include "liteseq/refs.h"	     // for ref_walk
#include "liteseq/types.h"	     // for strand
#include "povu/common/compat.hpp"    // for contains, pv_cmp
#include "povu/common/core.hpp"	     // for pt, idx_t, id_t, op_t
#include "povu/common/log.hpp"	     // for ERR
#include "povu/genomics/rov.hpp"     // for RoV
#include "povu/graph/bidirected.hpp" // for bd, VG
#include "povu/graph/pvst.hpp"	     // for VertexBase
#include "povu/graph/types.hpp"	     // for or_e, id_or_t, walk_t

namespace povu::genomics::allele
{
inline constexpr std::string_view MODULE = "povu::genomics::allele";

namespace lq = liteseq;
namespace pgt = povu::types::graph;
namespace pvst = povu::pvst;

struct allele_slice_t {
	const pgt::walk_t *walk;
	pt::idx_t walk_idx;
	pt::idx_t walk_start_idx;

	const lq::ref_walk *ref_w;
	pt::idx_t ref_idx;
	pt::idx_t ref_start_idx;

	pt::idx_t len; // total step count in the itinerary
	ptg::or_e slice_or;
	pgr::var_type_e vt;

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
		// auto [v_id, o_, _] = this->ref_walk->at(idx);
		return {ref_v_id, ref_o};
	}

	[[nodiscard]]
	std::string as_str() const
	{
		pt::idx_t ref_step_idx =
			this->get_or() == pgt::or_e::forward
				? this->ref_start_idx
				: this->ref_start_idx - len + 1;
		pt::idx_t end = ref_step_idx + this->len;
		std::string s;

		pt::u32 i = ref_step_idx;
		pt::u32 N = end;

		// INFO("allele_slice_t::as_str()");
		// std::cerr << " or " << this->get_or() << "\n";
		// std::cerr << " s " << this->ref_start_idx << " r idx "
		//	  << this->ref_idx << "\n";
		// std::cerr << "i " << i << " N " << N << "\n";

		switch (this->vt) {
		case pgr::var_type_e::sub:
			i++;
			N--;
			break;
		case pgr::var_type_e::ins:
		case pgr::var_type_e::del:
			N--;
			break;
		}

		for (i; i < N; i++)
			s += this->get_step(i).as_str();

		return s;
	}
};

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
	const pgr::RoV *rov_;

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

	Exp(const pgr::RoV *rov)
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
	const pgr::RoV *get_rov() const
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

std::vector<Exp> comp_itineraries3(const bd::VG &g, const pgr::RoV &rov,
				   const std::set<pt::id_t> &to_call_ref_ids);
std::vector<Exp> comp_itineraries2(const bd::VG &g, const pgr::RoV &rov);
void comp_itineraries(const bd::VG &g, Exp &exp);

} // namespace povu::genomics::allele

#endif // POVU_GENOMICS_ALLELE_HPP
