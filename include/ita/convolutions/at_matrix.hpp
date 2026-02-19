#include <thread>

#include <liteseq/refs.h> // for ref_walk, ref
#include <vector>

#include "ita/variation/rov.hpp"  // for RoV
#include "meza/matrix/matrix.hpp" // for matrix2d
#include "povu/common/core.hpp"
#include "povu/graph/bidirected.hpp" // for VG
#include "povu/graph/types.hpp"	     // for ptg: or_e

namespace ita::at_matrix
{

struct matrix_pool {
	std::map<pt::u32, meza::matrix::depth_matrix> ref_matrices;
	meza::matrix::depth_matrix filter_matrix;
	meza::matrix::depth_matrix result_matrix;
};

matrix_pool init_depth_matrices(const bd::VG &g, ir::RoV &rov,
				const std::set<pt::u32> &to_call_ref_ids);

} // namespace ita::at_matrix
