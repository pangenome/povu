#ifndef POVU_COLOR_HPP
#define POVU_COLOR_HPP

// #include <vector> // for vector
// #include "povu/common/app.hpp" // for config
//  #include "povu/common/constants.hpp"
#include "povu/common/core.hpp"	     // for pt constants
#include "povu/graph/bidirected.hpp" // for VG, bd
#include "povu/graph/pvst.hpp"	     // for Tree, VertexBase

// #include "povu/graph/types.hpp"	     // for or_e, id_or_t, walk_t

namespace povu::var::color
{
std::set<pt::u32> color_pvst(const bd::VG &g, const pvst::Tree &pvst,
			     const std::set<pt::id_t> &to_call_ref_ids);

} // namespace povu::var::color

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pvc = povu::var::color;

#endif // POVU_COLOR_HPP
