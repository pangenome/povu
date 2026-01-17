#ifndef ZIEN_HPP
#define ZIEN_HPP

#include <atomic> // for atomic
#include <ncurses.h>

#include "mto/from_vcf.hpp" // for VCFile

#include "povu/graph/bidirected.hpp" // for VG

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
		if (!isendwin()) {
			endwin();
		}
	}

	// Prevent copying to avoid multiple calls to endwin
	NcursesGuard(const NcursesGuard &) = delete;
	NcursesGuard &operator=(const NcursesGuard &) = delete;
};

void show_loading_spinner(std::atomic<bool> &is_loading);

void view(const bd::VG &g, const mto::from_vcf::VCFile &vcf_file,
	  const std::vector<pt::u32> &invalid_recs);
} // namespace zien::tui
#endif // ZIEN_HPP

namespace zt = zien::tui;
