#ifndef ITA_AT_MATRIX_HPP
#define ITA_AT_MATRIX_HPP

#include <optional>
#include <set>	  // for set
#include <vector> // for vector

#include <liteseq/refs.h>      // for ref_walk, ref
#include <meza/pool/split.hpp> // for matrix_pool

#include "ita/genomics/allele.hpp"	   // for hap_slice
#include "ita/traversals/depth_matrix.hpp" // for depth_matrix
#include "ita/traversals/traversals.hpp"   // for itinerary
#include "ita/variation/rov.hpp"	   // for RoV

#include "povu/common/core.hpp"	     // for pt
#include "povu/graph/bidirected.hpp" // for VG

#include "quilt/types.hpp"

namespace ita::at_matrix
{
using meza::pool::split::ov_mat_t;

struct hap2loop {
private:
	std::map<qt::u32, qt::u32> loop_pairing = {};

public:
	// ------------ constructors ----------

	hap2loop() = default;

	hap2loop(std::map<qt::u32, qt::u32> loop_pairing)
	    : loop_pairing(std::move(loop_pairing))
	{}

	// --------
	// accessor
	// --------

	[[nodiscard]]
	qt::u32 get_loop_no(qt::u32 h_idx) const
	{
		if (!qs::contains(loop_pairing, h_idx))
			return 0;

		return loop_pairing.at(h_idx);
	}
};

struct mat3 {
	ov_mat_t ref;
	ov_mat_t filter;
	ov_mat_t xor_result;

	qt::u32 j_offset = 0; // TODO: rename to pool j offset
	qt::u32 I = 0;
	qt::u32 J = 0;

	hap2loop h2l = {};

	// ------------
	// constructors
	// ------------

	mat3(ov_mat_t &&r, ov_mat_t &&f, ov_mat_t &&x, qt::u32 j_offset,
	     qt::u32 I, qt::u32 J)
	    : ref(std::move(r)), filter(std::move(f)), xor_result(std::move(x)),
	      j_offset(j_offset), I(I), J(J)
	{}

	mat3(ov_mat_t &&r, ov_mat_t &&f, ov_mat_t &&x, qt::u32 j_offset,
	     qt::u32 I, qt::u32 J, hap2loop &&h2l)
	    : ref(std::move(r)), filter(std::move(f)), xor_result(std::move(x)),
	      j_offset(j_offset), I(I), J(J), h2l(std::move(h2l))
	{}

	// disable copy constructor to prevent accidental copying
	mat3(const mat3 &other) = delete;
	mat3(mat3 &&other) noexcept = default;

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
	std::optional<tangle_info> tangle; // empty when no tangling

	// ------------
	// constructors
	// ------------

	mat3_item(mat3 &&m, std::optional<tangle_info> &&t)
	    : mats(std::move(m)), tangle(std::move(t))
	{}

	mat3_item(mat3 &&m) : mats(std::move(m)), tangle(std::nullopt)
	{}

	// ---------
	// accessors
	// ---------

	mat3 &get_mat3_mut()
	{
		return mats;
	}

	[[nodiscard]] bool is_tangled() const
	{
		return this->tangle != std::nullopt;
	}
};

struct rov_job {
private:
	using itn = ita::traversals::traversals::itinerary;

	const ir::RoV *rov = nullptr;

	// idx is ref_h_idx and value at idx is ...
	std::vector<std::vector<mat3_item>> items;

	// 1 is for no tangle, >1 is for tangled
	// the key is the ref idx of the mat3 item(s) ref
	// for any ref if the vector size is 1 there is no tangle,
	// >1 is for tangled
	// std::map<qt::u32, std::vector<mat3_item>> items;
	// std::vector<mat3_item> items;
	std::vector<itn> hap_itns;

public:
	// ------------
	// constructors
	// ------------

	rov_job() = delete;

	rov_job(const ir::RoV *r, qt::u32 H, std::vector<itn> hap_itns)
	    : rov(r), items(H), hap_itns(std::move(hap_itns))
	{}

	// --------
	// accessor
	// --------

	[[nodiscard]] const ir::RoV *get_rov() const
	{
		return rov;
	}

	[[nodiscard]] const std::vector<itn> &get_hap_itns() const
	{
		return hap_itns;
	}

	[[nodiscard]]
	const std::vector<mat3_item> &
	get_items_for_ref2(qt::u32 ref_h_idx) const
	{
		return items.at(ref_h_idx);
	}

	// --------
	// modifier
	// --------

	[[nodiscard]] std::vector<std::vector<mat3_item>> &get_items_mut()
	{
		return items;
	}

	void add_item(qt::u32 ref_h_idx, mat3_item &&item)
	{
		items[ref_h_idx].emplace_back(std::move(item));
	}
};

struct rov_job_batch {
private:
	std::vector<rov_job> items;
	qt::u32 pool_j_offset = 0;

	void comp_meta(qt::u32 j_idx)
	{
		rov_job &j = items.at(j_idx);
		for (std::vector<mat3_item> &ref_items : j.get_items_mut()) {
			for (mat3_item &ref_item_set : ref_items) {
				mat3 &m = ref_item_set.get_mat3_mut();
				m.ref.comp_meta_sync();
				m.filter.comp_meta_sync();
			}
		}
	}

public:
	// ------------
	// constructors
	// ------------

	// print when a copy or move is made for this struct
	rov_job_batch() = default;

	// disable copy constructor to prevent accidental copying
	rov_job_batch(const rov_job_batch &other) = delete;

	rov_job_batch(rov_job_batch &&other) noexcept
	    : items(std::move(other.items)), pool_j_offset(other.pool_j_offset)
	{}

	// ---------
	// accessors
	// ---------
	[[nodiscard]] qt::u32 get_pool_j_offset() const
	{
		return pool_j_offset;
	}

	[[nodiscard]] const std::vector<rov_job> &get_jobs()
	{
		return this->items;
	}

	// ---------
	// modifiers
	// ---------

	void reset()
	{
		items.clear();
		pool_j_offset = 0;
	}

	void set_pool_j_offset(qt::u32 offset)
	{
		pool_j_offset = offset;
	}

	void add(rov_job &&w)
	{
		items.emplace_back(std::move(w));
		std::thread([&]() { comp_meta(items.size() - 1); }).detach();
	}
};

void init_pool(const bd::VG &g, const ir::RoV *rov,
	       const std::set<pt::u32> &to_call_ref_ids,
	       const ita::depth_matrix::depth_matrix &dm,
	       meza::pool::split::matrix_pool<qt::u8> &ov_pool,
	       rov_job_batch &batch);
} // namespace ita::at_matrix
#endif // ITA_AT_MATRIX_HPP
