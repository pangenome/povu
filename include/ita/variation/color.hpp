#ifndef IT_COLOR_HPP
#define IT_COLOR_HPP

#include "povu/common/core.hpp"	     // for pt constants
#include "povu/graph/bidirected.hpp" // for VG, bd
#include "povu/graph/pvst.hpp"	     // for Tree, VertexBase

namespace ita::color
{
std::set<pt::u32> color_pvst(const bd::VG &g, const pvst::Tree &pvst,
			     const std::set<pt::id_t> &to_call_ref_ids);

} // namespace ita::color

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace ic = ita::color;

#endif // IT_COLOR_HPP
