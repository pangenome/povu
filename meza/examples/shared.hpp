#pragma once

#include "meza/pool/split.hpp"

namespace shared
{

void fill_random(meza::pool::split::ov_mat_t &mat);

void fill_ref_row(const meza::pool::split::ov_mat_t &filter_mat,
		  meza::pool::split::ov_mat_t &ref_mat);
} // namespace shared
