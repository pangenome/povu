#ifndef MT_IO_FRM_GFA_HPP
#define MT_IO_FRM_GFA_HPP

#include <string_view> // for string_view

#include <oza/graph/bidirected.hpp> // for VG
#include <quilt/app.hpp>	    // for config

namespace mto::from_gfa
{

/**
 * Remember to free the returned graph after use
 */
bd::VG *to_bd(const core::config &app_config);
}; // namespace mto::from_gfa
#endif // MT_IO_FRM_GFA_HPP
