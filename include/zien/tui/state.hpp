#ifndef ZIEN_STATE_HPP
#define ZIEN_STATE_HPP

#include <map>
#include <string>
#include <vector>

#include <ncurses.h>

#include "povu/common/core.hpp" // for pt

// Forward declaration: Tell the compiler Pane exists elsewhere

namespace zien::tui::state
{

// Constants for color pair IDs
static constexpr short PAIR_DEFAULT = 1;
static constexpr short PAIR_DIM = 2;
static constexpr short PAIR_SEP = 3;
static constexpr short PAIR_ERROR = 4;

enum class Mode : pt::u8 {
	NAVIGATION,
	SEARCH,
	JUMP,
	COMMAND,
};

enum class View : pt::u8 {
	VARIATION,
	PATHS,
};

extern std::map<Mode, std::string> mode_names;

enum class PaneID : uint8_t {
	// variation window
	A, // top left
	B, // top right
	C, // bottom left
	D, // bottom right
	E, // repeats pane
	// graph window
	F, // paths
};

extern std::map<PaneID, std::string> pane_names;

struct ui_state {

	// --------------------------
	// member variables or fields
	// --------------------------

	int screen_h; // screen height
	int screen_w; // screen width

	View current_view = View::VARIATION;

	PaneID active_pane_id = PaneID::A;

	pt::u8 pane_count = 4; // number of panes currently displayed

	// currently selected VCF record index
	pt::u32 vcf_selected_rec = 0; // zero-based index

	std::string jump_query = "";
	std::string search_query = "";
	std::string command_prompt = "";
	std::vector<int> search_results;
	int current_result_idx = -1;
	Mode current_mode = Mode::NAVIGATION;

	bool toggle_repeats_pane = false;

	// TODO: add is tangled to UI state

	// --------------------------
	// constructors & factory fns
	// ---------------------------

	// default constructor
	ui_state() = default;

	static ui_state create_new()
	{
		ui_state state{};
		getmaxyx(stdscr, state.screen_h, state.screen_w);

		return state;
	}

	// ----------------
	// methods
	// ----------------

	void setup_colors()
	{
		// -------------------------------------------------------------
		// colours
		//
		// ansi colours lookup
		// https://gist.github.com/JBlond/2fea43a3049b38287e5e9cefc87b2124
		// -------------------------------------------------------------

		if (!has_colors())
			return; // Safety check for terminals without color

		start_color();
		// use_default_colors(); // Optional: allows transparency in
		// some terminals

		short dark_gray_idx = 234;
		short medium_gray_idx = 240;

		init_pair(PAIR_DEFAULT, COLOR_WHITE, COLOR_BLACK);
		init_pair(PAIR_DIM, medium_gray_idx, COLOR_BLACK);
		init_pair(PAIR_SEP, dark_gray_idx, COLOR_BLACK);
		init_pair(PAIR_ERROR, COLOR_RED, COLOR_BLACK);
	}
};

}; // namespace zien::tui::state
#endif // ZIEN_STATE_HPP
