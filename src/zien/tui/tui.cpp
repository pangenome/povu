#include "zien/tui/tui.hpp"

#include <atomic>
#include <string>
#include <vector>

#include <ncurses.h>

#include <liteseq/refs.h> // for ref_walk, ref

#include "mto/from_vcf.hpp"		  // for VCFile
#include "povu/common/core.hpp"		  // for pt
#include "povu/graph/bidirected.hpp"	  // for VG
#include "zien/components/alts.hpp"	  // for update_alts
#include "zien/components/components.hpp" // for status_bar
#include "zien/components/genotypes.hpp"  // for update_haps
#include "zien/components/refs.hpp"	  // for foo
#include "zien/components/repeats.hpp"	  // for update_repeats
#include "zien/components/vcf.hpp"	  // for comp_vcfs
#include "zien/tui/state.hpp"		  // for Mode

namespace zien::tui
{
namespace lq = liteseq;

using namespace zien::tui::state;
using namespace zien::components;

const pt::u32 REPEATS_PANE_COUNT = 5;
const pt::u32 NO_REPEATS_PANE_COUNT = 4;

const pt::u32 INITIAL_VCF_REC_IDX = 0; // skip header line

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

void setup_pane(const pane_params &pp)
{
	auto [x, y, width, height, p_id, p, hh, sc, sl] = pp;

	p.win = newwin(height, width, y, x);
	p.width = width;
	p.height = height;
	p.p_id = p_id;
	p.title = pane_names[p_id];
	p.has_header = hh;
	p.sep_cols = sc;
}

void cleanup_panes(Pane &a, Pane &b, Pane &c, Pane &d, Pane &e)
{
	if (a.win)
		delwin(a.win);
	if (b.win)
		delwin(b.win);
	if (c.win)
		delwin(c.win);
	if (d.win)
		delwin(d.win);
	if (e.win)
		delwin(e.win);
}

void setup_layout(const ui_state &state, Pane &a, Pane &b, Pane &c, Pane &d,
		  Pane *e)
{
	int screen_h = state.screen_h;
	int screen_w = state.screen_w;
	// If e (Repeats) is NOT null, mid_h is screen_h / 3
	int mid_h = (e == nullptr) ? screen_h / 2 : screen_h / 3;
	int mid_w = screen_w / 2;

	int gutter = 1; // gutter of 1 char between windows

	setup_pane({0, 0, mid_w - gutter, mid_h, PaneID::A, a, true, false});
	setup_pane(
		{mid_w, 0, screen_w - mid_w, mid_h, PaneID::B, b, false, true});

	// 2. Bottom Left: Half width, bottom half
	// int bottom_h = screen_h - mid_h;
	// We subtract 1 from bottom_h if we want to leave the very last
	// line for the search bar
	int actual_view_h = mid_h - 1;

	actual_view_h = mid_h - 1; // keep bottom panes same height as top panes

	setup_pane({0, mid_h, mid_w - gutter, actual_view_h, PaneID::C, c});

	// 3. Bottom Right: Remaining width, bottom half
	setup_pane(
		{mid_w, mid_h, screen_w - mid_w, actual_view_h, PaneID::D, d});

	if (e != nullptr)
		setup_pane({0, mid_h + actual_view_h, screen_w, mid_h + 1,
			    PaneID::E, *e});
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

void show_loading_spinner(std::atomic<bool> &is_loading)
{
	curs_set(0); // hide the cursor

	const std::vector<std::string> frames = {"|", "/", "-", "\\"};
	int frame_idx = 0;

	// Set ncurses to non-blocking mode (wait 100ms for input)
	timeout(100);

	while (is_loading.load()) {
		erase();
		box(stdscr, 0, 0);

		// Fetch dimensions inside the loop in case the user
		// resizes the terminal
		int h, w;
		getmaxyx(stdscr, h, w);

		// No attron() or COLOR_PAIR needed for default white on
		// black
		mvprintw(h / 2, (w / 2) - 10, "%s Loading data...",
			 frames[frame_idx].c_str());

		refresh();

		// Cycle animation
		frame_idx = (frame_idx + 1) % frames.size();

		// Check for user input (optional: allow ESC to cancel)

		// Optional: allow user to break the loop manually
		if (getch() == 27) // ESC
			break;
	}

	// --- CLEANUP STEP ---
	erase();   // Clears the internal buffer (removes the box & text)
	refresh(); // Pushes that empty buffer to the physical screen

	curs_set(1); // restore the cursor

	// Reset ncurses to blocking mode for the rest of the app
	timeout(-1);
}

void update_status_bar(const mto::from_vcf::VCFile &vcf_file, ui_state &state,
		       status_bar &sb)
{
	// left
	//
	//
	// middle
	// const std::string &mode = mode_names.at(current_mode);
	const std::string &pane_name = pane_names.at(state.active_pane_id);

	// right
	std::string vcf_rec_line = pv_cmp::format(
		"[{}/{}]", state.vcf_selected_rec + 1, vcf_file.record_count());

	sb.draw(state, "", pane_name, vcf_rec_line);
};

void handle_special_states(int ch, tui_context &tc, ui_state &state)
{
	Pane *ap = tc.get_pane(state.active_pane_id); // active pane
	if (state.current_mode == Mode::SEARCH)
		handle_search_input(ch, ap, state);
	else if (state.current_mode == Mode::JUMP)
		handle_jump_input(ch, ap, state);
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

void draw_all_panes(const bd::VG &g, const mto::from_vcf::VCFile &vcf_file,
		    Pane &top_left_pane, Pane &top_right_pane,
		    Pane &bottom_left_pane, Pane &bottom_right_pane,
		    Pane &repeats_pane, ui_state &state, status_bar &sb,
		    bool initial_draw)
{

	if (initial_draw) {
		// Start at first record (skip header)
		top_left_pane.selected_line = 1;

		// rest
		top_right_pane.selected_line = 1;
		bottom_left_pane.selected_line = 0;
		bottom_right_pane.selected_line = 0;
	}

	pt::u32 rec_idx = top_left_pane.selected_line - 1;
	bool t = (vcf_file.get_records().at(rec_idx).is_tangled());

	if (t) {
		zien::components::repeats::update_repeats(g, vcf_file, rec_idx,
							  repeats_pane.pd);
	}

	top_left_pane.draw(state);
	top_right_pane.draw(state);
	bottom_left_pane.draw(state);
	bottom_right_pane.draw(state);
	update_status_bar(vcf_file, state, sb);

	if (t)
		repeats_pane.draw(state);

	if (initial_draw)
		refresh();
}

void update_depenent_panes(const bd::VG &g,
			   const mto::from_vcf::VCFile &vcf_file,
			   Pane &top_right_pane, Pane &bottom_left_pane,
			   Pane &bottom_right_pane, ui_state &state)
{
	top_right_pane.pd.lines.clear();
	zien::components::genotypes::update_haps(
		vcf_file, state.vcf_selected_rec, top_right_pane.pd);

	bottom_left_pane.pd.lines.clear();
	zien::components::refs::update_refs(g, vcf_file, state.vcf_selected_rec,
					    bottom_left_pane.pd);

	bottom_right_pane.pd.lines.clear();
	zien::components::alts::update_alts(g, vcf_file, state.vcf_selected_rec,
					    bottom_right_pane.pd);
}

void view(const bd::VG &g, const mto::from_vcf::VCFile &vcf_file,
	  const std::vector<pt::u32> &invalid_recs)
{
	tui_context tc;
	ui_state &state = tc.get_state();

	Pane *a = (tc.get_pane(PaneID::A));
	Pane *b = (tc.get_pane(PaneID::B));
	Pane *c = (tc.get_pane(PaneID::C));
	Pane *d = (tc.get_pane(PaneID::D));
	Pane *e = (tc.get_pane(PaneID::E));

	setup_layout(state, *a, *b, *c, *d, nullptr);

	status_bar sb(stdscr, state.screen_h - 1, state.screen_w);

	state.active_pane_id = PaneID::A; // Initial focus

	/* top left pane */
	// called once at the start
	zien::components::vcf::comp_vcfs(g, vcf_file, invalid_recs, a->pd);

	// 2. Explicitly trigger the data population BEFORE the loop
	b->pd.lines.clear();
	zien::components::genotypes::update_haps(vcf_file, INITIAL_VCF_REC_IDX,
						 b->pd);
	c->pd.lines.clear();
	zien::components::refs::update_refs(g, vcf_file, INITIAL_VCF_REC_IDX,
					    c->pd);
	d->pd.lines.clear();
	zien::components::alts::update_alts(g, vcf_file, INITIAL_VCF_REC_IDX,
					    d->pd);

	draw_all_panes(g, vcf_file, *a, *b, *c, *d, *e, state, sb, true);

	// aka has repeats
	// does the current record cover a repetitive region?
	bool is_tangled{false};

	int ch{};
	while ((ch = getch()) != 'q') {
		handle_input(ch, tc, state);

		// Update dependent panes if necessary
		if (state.active_pane_id == PaneID::A)
			update_depenent_panes(g, vcf_file, *b, *c, *d, state);

		// 2. Tangle Logic (Check if we need to change layout)
		is_tangled = vcf_file.get_records()
				     .at(state.vcf_selected_rec)
				     .is_tangled();

		state.pane_count =
			is_tangled ? REPEATS_PANE_COUNT : NO_REPEATS_PANE_COUNT;

		if (state.toggle_repeats_pane != is_tangled) {
			setup_layout(state, *a, *b, *c, *d,
				     (is_tangled ? e : nullptr));
			state.toggle_repeats_pane = is_tangled;

			erase();
			refresh();
		}

		draw_all_panes(g, vcf_file, *a, *b, *c, *d, *e, state, sb,
			       false);
		refresh(); // Push all changes to the physical terminal
	};

	return;
}

} // namespace zien::tui
