#ifndef ITA_AT_MATRIX_HPP
#define ITA_AT_MATRIX_HPP

#include <cstddef> // for size_t
#include <optional>
#include <set>	  // for set
#include <vector> // for vector

#include <convo/matrix.hpp> // for ref_matrix, depth_matrix, at_matrix
#include <convo/pool.hpp>   // for matrix_pool
#include <liteseq/refs.h>   // for ref_walk, ref

#include "ita/convolutions/depth_matrix.hpp" // for depth_matrix
#include "ita/genomics/allele.hpp"	     // for hap_slice
#include "ita/traversals/traversals.hpp"     // for itinerary
#include "ita/variation/rov.hpp"	     // for RoV
#include "povu/common/core.hpp"		     // for pt
#include "povu/graph/bidirected.hpp"	     // for VG
#include "quilt/types.hpp"

namespace ita::at_matrix
{

using ov_mat_t = meza::matrix_view::ov_matrix<qt::u8, std::string, std::string>;

struct mat3 {
	ov_mat_t ref;
	ov_mat_t filter;
	ov_mat_t xor_result;

	qt::u32 j_offset = 0; // TODO: rename to pool j offset
	qt::u32 I = 0;
	qt::u32 J = 0;

	void dbg_print() const
	{
		std::cerr << "Reference Matrix:\n";
		ref.dbg_print();
		std::cerr << "Filter Matrix:\n";
		filter.dbg_print();
		std::cerr << "XOR Result Matrix:\n";
		xor_result.dbg_print();
	}
};

struct tangle_info {
	qt::u32 ref_loop_no = 0;
	std::vector<std::optional<ia::hap_slice>> hap_slices;
};

struct mat3_item {
	mat3 mats;
	std::optional<tangle_info> tangle; // empty when untangled

	// ------------
	// constructors
	// ------------

	mat3_item(mat3 &&m, std::optional<tangle_info> &&t)
	    : mats(std::move(m)), tangle(std::move(t))
	{}

	mat3_item(mat3 &&m) : mats(std::move(m)), tangle(std::nullopt)
	{}
};

struct rov_job {
	const ir::RoV *rov = nullptr;
	// qt::u32 pool_j_offset = 0;

	// -------------
	// constructors
	// -------------

	// 1 is for no tangle, >1 is for tangled
	// the key is the ref idx of the mat3 item(s) ref
	// for any ref if the vector size is 1 there is no tangle,
	// >1 is for tangled
	std::map<qt::u32, std::vector<mat3_item>> items;
	// std::vector<mat3_item> items;
	std::vector<ita::traversals::traversals::itinerary> hap_itns;

	rov_job() = delete;

	rov_job(const ir::RoV *r,
		std::vector<ita::traversals::traversals::itinerary> hap_itns)
	    : rov(r), hap_itns(std::move(hap_itns))
	{}

	void add_item(qt::u32 ref_h_idx, mat3_item &&item)
	{
		items[ref_h_idx].emplace_back(std::move(item));
	}

	// rov_job(const ir::RoV *r, qt::u32 j_offset, qt::u32 ref_h_idx,
	//	mat3_item &&item)
	//     : rov(r), pool_j_offset(j_offset), hap_itns({})
	// {
	//	items[ref_h_idx].emplace_back(std::move(item));
	// }

	// rov_job(const ir::RoV *r, qt::u32 j_offset, qt::u32 ref_h_idx,
	//	const std::vector<mat3_item> &items,
	//	std::vector<ita::traversals::traversals::itinerary> hap_itns)
	//     : rov(r), pool_j_offset(j_offset), hap_itns(std::move(hap_itns))
	// {
	//	std::vector<mat3_item> &ref_items = this->items[ref_h_idx];
	//	ref_items.insert(ref_items.end(), items.begin(), items.end());
	// }

	// -------
	// setters
	// -------

	[[nodiscard]]
	const std::vector<mat3_item> *get_items_for_ref(qt::u32 ref_h_idx) const
	{
		auto it = items.find(ref_h_idx);

		if (it == items.end())
			return nullptr;

		return &it->second;
	}
};

struct rov_job_batch {
	std::vector<rov_job> items;
	qt::u32 pool_j_offset = 0;

	// [[nodiscard]]
	// bool empty() const
	// {
	//	return items.empty();
	// }

	// [[nodiscard]]
	// std::size_t size() const
	// {
	//	return items.size();
	// }

	// rov_job &back()
	// {
	//	return items.back();
	// }

	// [[nodiscard]]
	// const rov_job &at(std::size_t i) const
	// {
	//	return items[i];
	// }

	void reset()
	{
		items.clear();
		pool_j_offset = 0;
	}

	void add(rov_job &&w)
	{
		items.emplace_back(std::move(w));
	}
};

void init_pool(const bd::VG &g, const ir::RoV *rov,
	       const std::set<pt::u32> &to_call_ref_ids,
	       const ita::depth_matrix::depth_matrix &dm,
	       meza::matrix_pool::matrix_pool<qt::u8> &ov_pool,
	       rov_job_batch &batch);

// struct rov_mat_set {
//	ov_mat_t ref;
//	ov_mat_t filter;
//	ov_mat_t xor_result;

//	qt::u32 mat_j_offset;
//	const ir::RoV *rov;

//	void dbg_print() const
//	{
//		std::cerr << "Reference Matrix:\n";
//		ref.dbg_print();
//		std::cerr << "Filter Matrix:\n";
//		filter.dbg_print();
//		std::cerr << "XOR Result Matrix:\n";
//		xor_result.dbg_print();
//	}
// };

// struct all_mat_sets {
//	std::vector<rov_mat_set> mat_sets;
//	qt::u32 pool_j_offset = 0;

//	// -------
//	// getters
//	// -------

//	[[nodiscard]]
//	bool empty() const
//	{
//		return mat_sets.empty();
//	}

//	[[nodiscard]]
//	std::size_t size() const
//	{
//		return this->mat_sets.size();
//	}

//	[[nodiscard]]
//	const ita::at_matrix::rov_mat_set &at(qt::u32 i) const
//	{
//		return mat_sets[i];
//	}

//	// ---------
//	// modifiers
//	// ---------

//	void reset()
//	{
//		this->mat_sets.clear();
//		this->pool_j_offset = 0;
//	}

//	void add(ita::at_matrix::rov_mat_set &&s)
//	{
//		mat_sets.emplace_back(std::move(s));
//	}
// };

// struct mat_bundle {

//	const ir::RoV *rov;

//	ov_mat_t ref;
//	ov_mat_t filter;
//	ov_mat_t xor_result;
//	qt::u32 mat_j_offset;

//	// only used for tangled pools, otherwise 0 by default
//	pt::u32 ref_loop_no = 0;
//	std::vector<std::optional<ia::hap_slice>> hap_slices;

//	void dbg_print() const
//	{
//		std::cerr << "Reference Matrix:\n";
//		ref.dbg_print();
//		std::cerr << "Filter Matrix:\n";
//		filter.dbg_print();
//		std::cerr << "XOR Result Matrix:\n";
//		xor_result.dbg_print();
//	}
// };

// struct mat_triple {
//	ov_mat_t ref;
//	ov_mat_t filter;
//	ov_mat_t xor_result;
//	qt::u32 mat_j_offset;
//	pt::u32 I;
//	pt::u32 J;

//	mat_triple(ov_mat_t &&r, ov_mat_t &&f, ov_mat_t &&x,
//		   qt::u32 mat_j_offset, pt::u32 I, pt::u32 J)
//	    : ref(std::move(r)), filter(std::move(f)), xor_result(std::move(x)),
//	      mat_j_offset(mat_j_offset), I(I), J(J)
//	{}
// };

// struct mat_triple_tangled : public mat_triple {
//	pt::u32 ref_loop_no = 0;
//	std::vector<std::optional<ia::hap_slice>> hap_slices;

//	mat_triple_tangled(
//		ov_mat_t &&r, ov_mat_t &&f, ov_mat_t &&x, qt::u32 mat_j_offset,
//		pt::u32 I, pt::u32 J, pt::u32 ref_loop_no,
//		std::vector<std::optional<ia::hap_slice>> &&hap_slices)
//	    : mat_triple(std::move(r), std::move(f), std::move(x), mat_j_offset,
//			 I, J),
//	      ref_loop_no(ref_loop_no), hap_slices(std::move(hap_slices))
//	{}
// };

// struct bundle {
//	bool is_tangled_;
//	const ir::RoV *rov;
//	qt::u32 pool_j_offset;
//	std::vector<ita::traversals::traversals::itinerary> hap_itns_;
//	std::variant<mat_triple, std::vector<mat_triple_tangled>> mat_data_;

//	bundle(bool is_tangled, const ir::RoV *r, qt::u32 j_offset,
//	       std::vector<ita::traversals::traversals::itinerary> hap_itns,
//	       std::vector<mat_triple_tangled> &&tangled_data)
//	    : is_tangled_(is_tangled), rov(r), pool_j_offset(j_offset),
//	      hap_itns_(std::move(hap_itns)), mat_data_(std::move(tangled_data))
//	{}

//	[[nodiscard]] bool is_tangled() const
//	{
//		return is_tangled_;
//	}

//	std::vector<mat_triple_tangled> &get_tangled_mats()
//	{
//		return std::get<std::vector<mat_triple_tangled>>(mat_data_);
//	}

//	void set_tangled(bool t)

//	{
//		is_tangled_ = t;
//	}
// };

// struct rov_mat_bundles {
//	bool tangled = false; // if not tangled bundles should have size 1
//	std::vector<mat_bundle> bundles;
//	std::vector<ita::traversals::traversals::itinerary> hap_itns;

//	rov_mat_bundles(
//		bool t,
//		const std::vector<ita::traversals::traversals::itinerary>
//			&hap_itns)
//	    : tangled(t), hap_itns(hap_itns)
//	{}
// };

// struct pool_mat_bundles {
//	std::vector<rov_mat_bundles> data_;
//	qt::u32 pool_j_offset = 0;

//	void add(rov_mat_bundles &&e)
//	{
//		data_.emplace_back(std::move(e));
//	}

//	rov_mat_bundles &back()
//	{
//		return data_.back();
//	}
// };

// struct matrix_pool {

//	std::map<pt::u32, meza::matrix::ref_matrix> ref_matrices;
//	meza::matrix::depth_matrix filter_matrix;
//	meza::matrix::depth_matrix result_matrix;
//	std::vector<std::optional<ia::hap_slice>> hap_slices;

//	// only used for tangled pools, otherwise 0 by default
//	pt::u32 ref_loop_no = 0;

//	bool tangled;
//	pt::u32 I;
//	pt::u32 J;
// };

// struct rov_matrix_pool {
//	const ir::RoV &rov;
//	std::vector<matrix_pool> pools;
//	bool tangled = false; // if not tangled pools should have size 1
//	std::vector<ita::traversals::traversals::itinerary> hap_itns;

//	// --------------
//	// constructor(s)
//	// --------------

//	rov_matrix_pool(const ir::RoV &rov) : rov(rov)
//	{}

//	rov_matrix_pool(
//		const ir::RoV &rov, bool t,
//		const std::vector<ita::traversals::traversals::itinerary>
//			&hap_itns)
//	    : rov(rov), tangled(t), hap_itns(hap_itns)
//	{}

//	// ---------
//	// getter(s)
//	// ---------

//	[[nodiscard]]
//	bool empty() const
//	{
//		return pools.empty();
//	}
// };

// rov_matrix_pool init_depth_matrices(const bd::VG &g, const ir::RoV *rov,
//				    const std::set<pt::u32> &to_call_ref_ids);
} // namespace ita::at_matrix
#endif // ITA_AT_MATRIX_HPP
