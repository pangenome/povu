#ifndef ZIEN_COMPONENTS_PATHS_HPP
#define ZIEN_COMPONENTS_PATHS_HPP

#include <ncurses.h>

#include "povu/graph/bidirected.hpp" // for VG

namespace zien::components::paths
{
void update_paths(const bd::VG &g, display_lines &pd);
} // namespace zien::components::paths

#endif // ZIEN_COMPONENTS_PATHS_HPP
