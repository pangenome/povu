#pragma once

#include <liteseq/refs.h>	 // for ref_walk, ref
#include <meza/owned/matrix.hpp> // for dense_matrix2d
#include <meza/pool/joint.hpp>	 // for joint_pool
#include <meza/pool/split.hpp>	 // for matrix_pool

#include "ita/variation/rov.hpp"     // for RoV
#include "povu/graph/bidirected.hpp" // for VG

namespace ita::depth_matrix
{

struct depth_matrix {
private:
	bool tangled_ = false;
	qt::u32 max_depth_ = 0;
	meza::pool::joint::full_view<qt::u32> view;

public:
	// -----------
	// constructor
	// -----------

	depth_matrix(meza::pool::joint::joint_pool<qt::u32> &pool, qt::u32 I,
		     qt::u32 J)
	    : view(pool.alloc_full(I, J))
	{}

	// -------
	// getters
	// -------

	[[nodiscard]]
	qt::u32 rows() const
	{
		return view.rows();
	}

	[[nodiscard]]
	qt::u32 cols() const
	{
		return view.cols();
	}

	[[nodiscard]]
	std::pair<const qt::u32 *, qt::u32> get_slice() const
	{
		return this->view.get_slice();
	}

	[[nodiscard]]
	const qt::u32 &get_value(qt::u32 i, qt::u32 j) const
	{
		return this->view.at(i, j);
	}

	[[nodiscard]]
	const qt::u32 &get_max_depth() const
	{
		return this->max_depth_;
	}

	[[nodiscard]]
	bool is_tangled() const
	{
		return tangled_;
	}

	// -------
	// setters
	// -------

	void set_max_depth(qt::u32 max_depth)
	{
		this->max_depth_ = max_depth;
	}

	void set_value(qt::u32 i, qt::u32 j, qt::u32 val)
	{
		this->view.at_mut(i, j) = val;
	}

	void set_tangled(bool tangled)
	{
		tangled_ = tangled;
	}
};

depth_matrix comp_depth_matrix(const bd::VG &g, const ir::RoV *rov,
			       meza::pool::joint::joint_pool<qt::u32> &dm_pool);

}; // namespace ita::depth_matrix
