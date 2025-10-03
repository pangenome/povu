#ifndef PV_IO_FRM_GFA_HPP
#define PV_IO_FRM_GFA_HPP

"#include "povu/common/app.hpp"
"#include "povu/graph/bidirected.hpp"

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
