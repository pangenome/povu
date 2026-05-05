#ifndef IT_COLOR_HPP
#define IT_COLOR_HPP

#include <oza/graph/bidirected.hpp> // for VG, bd
#include <oza/graph/pvst.hpp>	    // for Tree, VertexBase
#include <quilt/types.hpp>	    // for qt

namespace ita::color
{
std::set<qt::u32> color_pvst(const bd::VG &g, const pvst::Tree &pvst,
			     const std::set<qt::id_t> &to_call_ref_ids);

} // namespace ita::color

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace ic = ita::color;

#endif // IT_COLOR_HPP
