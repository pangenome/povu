#ifndef ZIEN_HPP
#define ZIEN_HPP

#include <atomic> // for atomic
#include <ncurses.h>

#include "mto/from_vcf.hpp"	     // for VCFile
#include "povu/graph/bidirected.hpp" // for VG
#include "zien/components/components.hpp"
#include "zien/tui/state.hpp" // for Mode

namespace zien::tui
{
inline constexpr std::string_view MODULE = "zien";

struct NcursesGuard {
	NcursesGuard()
	{
		initscr();	      // Initialize ncurses
		noecho();	      // Don't echo keypresses
		cbreak();	      // Disable line buffering
		keypad(stdscr, TRUE); // Enable arrow keys/F-keys
		if (has_colors()) {
			start_color();
			use_default_colors(); // Allows COLOR_BLACK to be
					      // transparent if desired
		}
	}

	// The destructor runs automatically when the function returns or
	// crashes
	~NcursesGuard()
	{
		if (!isendwin())
			endwin();
	}

	// Prevent copying to avoid multiple calls to endwin
	NcursesGuard(const NcursesGuard &) = delete;
	NcursesGuard &operator=(const NcursesGuard &) = delete;
};

struct tui_context {
	// -------------
	// member fields
	// -------------
	zien::tui::state::ui_state state;

	// 1. Store the actual Pane objects
	components::Pane top_left_pane;
	components::Pane top_right_pane;
	components::Pane bottom_left_pane;
	components::Pane bottom_right_pane;
	components::Pane repeats_pane;
	components::Pane paths_pane;

	std::map<zien::tui::state::PaneID, components::Pane *> panes;

	// -----------
	// constructor
	// -----------
	tui_context(const tui_context &) = delete;	      // No copying
	tui_context &operator=(const tui_context &) = delete; // No assignment

	tui_context()
	    : state(zien::tui::state::ui_state::create_new()),
	      panes{{zien::tui::state::PaneID::A, &top_left_pane},
		    {zien::tui::state::PaneID::B, &top_right_pane},
		    {zien::tui::state::PaneID::C, &bottom_left_pane},
		    {zien::tui::state::PaneID::D, &bottom_right_pane},
		    {zien::tui::state::PaneID::E, &repeats_pane},
		    {zien::tui::state::PaneID::F, &paths_pane}}
	{
		state.setup_colors();
	}

	// -------
	// methods
	// -------

	zien::tui::state::ui_state &get_state()
	{
		return this->state;
	}

	components::Pane &get_pane_ref(zien::tui::state::PaneID pane_id)
	{
		return *this->panes.at(pane_id);
	}

	components::Pane *get_pane(zien::tui::state::PaneID pane_id)
	{
		return this->panes.at(pane_id);
	}
};

void show_loading_spinner(std::atomic<bool> &is_loading);

void view_gfa(const bd::VG &g);

void view(const bd::VG &g, const mto::from_vcf::VCFile &vcf_file,
	  const std::vector<pt::u32> &invalid_recs);
} // namespace zien::tui
#endif // ZIEN_HPP

namespace zt = zien::tui;
