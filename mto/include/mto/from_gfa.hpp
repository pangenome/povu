#ifndef MT_IO_FRM_GFA_HPP
#define MT_IO_FRM_GFA_HPP

#include <string_view> // for string_view

#include "povu/common/app.hpp"	     // for config
#include "povu/graph/bidirected.hpp" // for VG

namespace mto::from_gfa
{
inline constexpr std::string_view MODULE = "povu::io::from_gfa";
namespace bd = povu::bidirected;

/**
 * Remember to free the returned graph after use
 */
bd::VG *to_bd(const core::config &app_config);
}; // namespace mto::from_gfa
#endif // MT_IO_FRM_GFA_HPP
