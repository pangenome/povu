#pragma once

#include <liteseq/refs.h> // for ref_walk, ref
// #include <meza/owned/matrix.hpp> // for dense_matrix2d
#include <meza/pool/joint.hpp> // for joint_pool
// #include <meza/pool/split.hpp>	 // for matrix_pool

#include <meza/pool/pool.hpp> // for pool

#include "ita/variation/rov.hpp"     // for RoV
#include "povu/graph/bidirected.hpp" // for VG
#include "quilt/types.hpp"

namespace ita::depth_matrix
{
using pool_t = meza::pool::pool<qt::u8, qt::u32>;

struct depth_matrix {
private:
	bool tangled_ = false;
	qt::u32 max_depth_ = 0;
	meza::pool::joint::full_view<qt::u32> view;

public:
	// -----------
	// constructor
	// -----------

	depth_matrix(meza::pool::joint::full_view<qt::u32> &&v) : view(v)
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

depth_matrix comp_depth_matrix(const bd::VG &g, const ir::RoV *rov, pool_t &p);

}; // namespace ita::depth_matrix
