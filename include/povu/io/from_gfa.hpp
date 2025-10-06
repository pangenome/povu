#ifndef PV_IO_FRM_GFA_HPP
#define PV_IO_FRM_GFA_HPP

#include <string_view> // for string_view

#include "povu/common/app.hpp"	     // for config
#include "povu/graph/bidirected.hpp" // for VG

namespace povu::io::from_gfa
{
inline constexpr std::string_view MODULE = "povu::io::from_gfa";
namespace bd = povu::bidirected;

/**
 * Remember to free the returned graph after use
 */
bd::VG *to_bd(const core::config &app_config);
}; // namespace povu::io::from_gfa
#endif // PV_IO_FRM_GFA_HPP
