#include "zien/tui/state.hpp" // for Mode

#include <ncurses.h>

namespace zien::tui::state
{
std::map<Mode, std::string> mode_names = {
	{Mode::NAVIGATION, "NAVIGATION"},
	{Mode::SEARCH, "SEARCH"},
	{Mode::JUMP, "JUMP"},
};

std::map<PaneID, std::string> pane_names = {{PaneID::A, "VCF"},
					    {PaneID::B, "Haps"},
					    {PaneID::C, "Refs"},
					    {PaneID::D, "Alts"},
					    {PaneID::E, "Repeats"}};

}; // namespace zien::tui::state
