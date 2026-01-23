#include <string>
#include <vector>

#include <ncurses.h>

#include <liteseq/refs.h> // for ref_walk, ref

#include "povu/common/core.hpp"		  // for pt
#include "zien/components/components.hpp" // for status_bar
#include "zien/tui/state.hpp"		  // for Mode
#include "zien/tui/tui.hpp"

namespace zien::tui::input
{
using namespace zien::components;
using namespace zien::tui::state;

void perform_search(Pane *pane, ui_state &state)
{
	std::string &query = state.search_query;
	state.search_results.clear();
	if (query.empty())
		return;

	for (pt::u32 i = 0; i < pane->pd.lines.size(); ++i) {
		if (pane->pd.lines[i].find(query) != std::string::npos) {
			state.search_results.push_back(i);
		}
	}

	if (!state.search_results.empty()) {
		// Find the first result that is >= the currently
		// selected line
		state.current_result_idx = 0;
		for (size_t i = 0; i < state.search_results.size(); ++i) {
			if (state.search_results[i] >= pane->selected_line) {
				state.current_result_idx = i;
				break;
			}
		}
		pane->selected_line =
			state.search_results[state.current_result_idx];
		pane->horiz_offset = 0;
		pane->update_scroll();
	}
	else {
		state.current_result_idx = -1;
	}
}

void clear_search(ui_state &state)
{
	state.search_results.clear();
	state.current_result_idx = -1;
}

// handle search input
void handle_search_input(int ch, Pane *ap, ui_state &state)
{
	switch (ch) {
	case '\r':
	case '\n':
		state.current_mode = Mode::NAVIGATION;
		perform_search(ap, state);
		break;
	case 27: // ESC key
		state.current_mode = Mode::NAVIGATION;
		clear_search(state);
		break;
	case KEY_BACKSPACE:
	case 127: // Handle Backspace or DEL
		if (!state.search_query.empty())
			state.search_query.pop_back();
		break;
	default:
		if (ch >= 32 && ch <= 126) // Only add printable characters
			state.search_query += (char)ch;
	}
}

void handle_jump_input(int ch, Pane *active, ui_state &state)
{
	// Pane *active = ;
	switch (ch) {
	case '\n': // ENTER: Execute jump
		if (!state.jump_query.empty()) {
			try {
				int target = std::stoi(state.jump_query);
				int min_idx = active->has_header ? 1 : 0;
				int max_idx = (int)active->pd.lines.size() - 1;

				// Bounds checking
				if (target >= min_idx && target <= max_idx) {
					active->selected_line = target;
				}
				else if (target > max_idx) {
					active->selected_line = max_idx;
				}
				else {
					active->selected_line = min_idx;
				}
				active->update_scroll();
			}
			catch (...) {
			} // Handle non-numeric junk
		}
		state.current_mode = Mode::NAVIGATION;
		state.jump_query = "";
		break;

	case 27: // ESC: Cancel
		state.current_mode = Mode::NAVIGATION;
		state.jump_query = "";
		break;

	case KEY_BACKSPACE:
	case 127:
		if (!state.jump_query.empty())
			state.jump_query.pop_back();
		break;

	default:
		if (isdigit(ch)) { // Only allow numbers for jumping
			state.jump_query += (char)ch;
		}
		break;
	}
}

void commands_handler(ui_state &state)
{
	const std::string &cmd = state.command_prompt;
	if (cmd == "q" || cmd == "quit") {
		// Exit the application
		endwin(); // End ncurses mode
		std::exit(0);
	}
	else if (cmd == "h" || cmd == "help") {
		// TODO
	}
	else if (cmd == "g") {
		state.current_view = zien::tui::state::View::PATHS;
		state.active_pane_id = PaneID::F;
	}
	else if (cmd == "v") {
		state.current_view = zien::tui::state::View::VARIATION;
		state.active_pane_id = PaneID::A;
	}
}

// handle search input
void handle_command_input(int ch, Pane *ap, ui_state &state)
{
	switch (ch) {
	case '\r':
	case '\n':
		commands_handler(state);
		state.current_mode = Mode::NAVIGATION;
		state.command_prompt.clear();
		break;
	case 27: // ESC key
		state.current_mode = Mode::NAVIGATION;
		state.command_prompt.clear();
		break;
	case KEY_BACKSPACE:
	case 127: // Handle Backspace or DEL
		if (!state.command_prompt.empty())
			state.command_prompt.pop_back();
		break;
	default:
		if (ch >= 32 && ch <= 126) // Only add printable characters
			state.command_prompt += (char)ch;
	}
}

void handle_special_states(int ch, tui_context &tc, ui_state &state)
{
	Pane *ap = tc.get_pane(state.active_pane_id); // active pane
	if (state.current_mode == Mode::SEARCH)
		handle_search_input(ch, ap, state);
	else if (state.current_mode == Mode::JUMP)
		handle_jump_input(ch, ap, state);
	else if (state.current_mode == Mode::COMMAND)
		handle_command_input(ch, ap, state);
}

// navigate and initialize search
void nav(int ch, tui_context &tc, ui_state &state)
{
	PaneID active_pane_id = state.active_pane_id;
	Pane *active = tc.get_pane(active_pane_id);

	switch (ch) {
	case 27: // ESC key
		state.current_mode = Mode::NAVIGATION;
		clear_search(state);
		break;
	case '/':
		state.current_mode = Mode::SEARCH;
		state.search_query = "";
		break;
	case ':':
		state.current_mode = Mode::JUMP;
		break;
	case ' ':
		state.current_mode = Mode::COMMAND;
		break;
	case 'n': // Next match (Forward)
		if (!state.search_results.empty()) {
			state.current_result_idx =
				(state.current_result_idx + 1) %
				state.search_results.size();
			active->selected_line =
				state.search_results[state.current_result_idx];
		}
		break;
	case 'N': // Previous match (Backward)
		if (!state.search_results.empty()) {
			// Adding size before modulo handles the
			// negative wrap-around
			state.current_result_idx =
				(state.current_result_idx - 1 +
				 state.search_results.size()) %
				state.search_results.size();
			active->selected_line =
				state.search_results[state.current_result_idx];
			active->update_scroll();
		}
		break;
	case '\t': // Cycle: Top -> Left -> Right -> Top
		switch (active_pane_id) {
		case PaneID::A:
			state.active_pane_id = PaneID::B;
			break;
		case PaneID::B:
			state.active_pane_id = PaneID::C;
			break;
		case PaneID::C:
			state.active_pane_id = PaneID::D;
			break;
		case PaneID::D:
			state.active_pane_id =
				state.pane_count == 5 ? PaneID::E : PaneID::A;
			break;
		case PaneID::E:
			state.active_pane_id = PaneID::A;
			break;
		default:
			state.active_pane_id = PaneID::A;
		}
		break;
	case 'j':
	case KEY_DOWN:
		active->scroll_down();
		break;
	case 'k':
	case KEY_UP:
		active->scroll_up();
		break;
	case 'l':
	case KEY_RIGHT:
		if (active->horiz_offset + active->get_content_width() <
		    (int)active->pd.lines[active->selected_line].length())
			active->horiz_offset++;
		break;
	case 'h':
	case KEY_LEFT:
		if (active->horiz_offset > 0)
			active->horiz_offset--;
		break;
	}

	active->update_scroll();
}

void handle_normal_state(int ch, tui_context &tc, ui_state &state)
{
	switch (ch) {
	case '/':
		state.current_mode = Mode::SEARCH;
		break;
	case ':':
		state.current_mode = Mode::JUMP;
		break; // Use colon for jump
	case ' ':
		// Toggle repeats pane visibility
		state.current_mode = Mode::COMMAND;
		break;
	}

	// nav & search
	nav(ch, tc, state);

	// Immediately sync state ONLY if the VCF pane was the
	// one moving
	if (state.active_pane_id == PaneID::A) {
		Pane *a = tc.get_pane(PaneID::A);
		state.vcf_selected_rec = (pt::u32)a->selected_line - 1;
	}
}

void handle_input(int ch, tui_context &tc, ui_state &state)
{
	if (state.current_mode == Mode::NAVIGATION)
		handle_normal_state(ch, tc, state);
	else
		handle_special_states(ch, tc, state);
}
}; // namespace zien::tui::input
