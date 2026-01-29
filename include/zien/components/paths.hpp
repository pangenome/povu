#ifndef ZIEN_COMPONENTS_PATHS_HPP
#define ZIEN_COMPONENTS_PATHS_HPP

#include <ncurses.h>

#include "povu/graph/bidirected.hpp"	  // for VG
#include "zien/components/components.hpp" // for display_lines
#include "zien/tui/state.hpp"		  // for ui_state

namespace zien::components::paths
{
void update_paths(const bd::VG &g, ui_state &state, display_lines &pd);
} // namespace zien::components::paths

#endif // ZIEN_COMPONENTS_PATHS_HPP
