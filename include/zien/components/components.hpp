#ifndef ZIEN_COMPONENTS_HPP
#define ZIEN_COMPONENTS_HPP

#include <iostream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

#include <ncurses.h>

#include "povu/common/compat.hpp" // for pv_cmp
#include "zien/tui/state.hpp"	  // for Mode

namespace zien::components
{
using namespace zien::tui::state;

struct line_metadata {
	pt::u32 ref_name_pos;
	std::vector<pt::slice> at_str_slices; // allele traversal string slices
};

// the lines to display in a given pane
struct display_lines {
	std::vector<std::string> lines;	 // all lines to display
	std::set<pt::u32> special_lines; // lines to highlight (error lines)
	std::map<pt::u32, line_metadata> meta; // all lines metadata
	std::set<pt::u32> group_lines;	       // lines that start new groups
	pt::u32 lh = 0;			       // longest header or label width

	void reset()
	{
		this->lines.clear();
		this->meta.clear();
		this->group_lines.clear();
		this->lh = 0;
	}
};

struct Pane {
	WINDOW *win;
	std::string title;	 // title field
	bool has_header = false; // only set true for top pane, for now at least

	int width;
	int height;

	PaneID p_id;	       // pane identifier
	bool sep_cols = false; // whether to separate columns (for tabular data)

	int scroll_offset = 0;
	int horiz_offset = 0; // Horizontal scroll position
	int selected_line = 0;

	display_lines pd;

	// -------------------------
	// constructors & destructor
	// -------------------------
	Pane() = default;
	Pane(const Pane &) = delete;		// No copying allowed
	Pane &operator=(const Pane &) = delete; // No assignment allowed

	~Pane()
	{
		if (win)
			delwin(win);
	}

	// ----------
	// methods
	// ----------
	int get_content_width()
	{
		int safe_w = width - 2;
		if (p_id == PaneID::C || p_id == PaneID::D) {
			int label_width = pd.lh;
			if (label_width > safe_w * 0.4)
				label_width = safe_w * 0.4;
			return safe_w - label_width - 2; // 2 for gutter
		}
		return safe_w;
	}

	void update_scroll()
	{
		int data_area_h = height - 2 - (has_header ? 2 : 0);
		int first_data_idx = has_header ? 1 : 0;

		// Calculate how many EXTRA visual rows exist due to separators
		// from the start of the data up to the selected line.
		int separators_before_selection = 0;
		for (auto line_idx : pd.group_lines) {
			// Only count separators that are within the current
			// scroll view OR affect the current selection's
			// position.
			if ((int)line_idx >= first_data_idx &&
			    (int)line_idx <= selected_line) {
				separators_before_selection++;
			}
		}

		// This is the actual screen row (relative to data start)
		// the selection wants to occupy.
		int virtual_selected_pos = (selected_line - first_data_idx) +
					   separators_before_selection;

		if (virtual_selected_pos < scroll_offset) {
			scroll_offset = virtual_selected_pos;
		}

		if (virtual_selected_pos >= scroll_offset + data_area_h) {
			scroll_offset = virtual_selected_pos - data_area_h + 1;
		}

		if (scroll_offset < 0)
			scroll_offset = 0;
	}

	void draw_tabular_line(int y, int x, const std::string &line,
			       int col_width)
	{
		std::stringstream ss(line);
		std::string segment;
		int current_col = 0;

		// Split the string by tabs
		while (std::getline(ss, segment, '\t')) {
			int start_x = x + (current_col * col_width);

			// Truncate segment if it's wider than the column to
			// prevent overlap
			if (segment.length() > (size_t)col_width - 1) {
				segment =
					segment.substr(0, col_width - 2) + "~";
			}

			mvwprintw(win, y, start_x, "%s", segment.c_str());
			current_col++;
		}
	}

	void scroll_down()
	{
		if (pd.lines.empty())
			return;

		pt::u32 min_idx = has_header ? 1 : 0;
		pt::u32 max_idx = pd.lines.size() - 1;

		if ((pt::u32)selected_line < max_idx)
			selected_line++;
		else
			selected_line = min_idx;
	}

	void scroll_up()
	{
		if (pd.lines.empty())
			return;

		pt::u32 min_idx = has_header ? 1 : 0;
		pt::u32 max_idx = pd.lines.size() - 1;

		if ((pt::u32)selected_line > min_idx)
			selected_line--;
		else
			selected_line = max_idx;
	}

	void draw(const ui_state &state)
	{
		bool is_focused = (state.active_pane_id == this->p_id);

		update_scroll(); // Ensure this uses the "visual row" logic
				 // discussed
		werase(win);
		box(win, 0, 0);

		const std::vector<std::string> &lines = this->pd.lines;
		if (lines.empty()) {
			box(win, 0, 0);
			wrefresh(win);
			return;
		}

		if (!title.empty())
			mvwprintw(win, 0, 2, " %s ", title.c_str());

		int safe_w = width - 2;
		int data_start_y = has_header ? 3 : 1;
		int body_h = (height - 2) - (has_header ? 2 : 0);
		int first_data_idx = has_header ? 1 : 0;

		// --- 1. HEADER RENDERING ---
		if (has_header) {
			wattron(win, A_BOLD | COLOR_PAIR(2));
			mvwprintw(win, 1, 1, "%-*.*s", safe_w, safe_w,
				  lines[0].c_str());
			wattroff(win, A_BOLD | COLOR_PAIR(2));
			mvwhline(win, 2, 1, ACS_HLINE, safe_w);
		}

		// --- 2. COLUMN & LABEL WIDTH CALCULATIONS ---
		int col_width = 10;
		if (sep_cols) {
			// Sample visible lines to determine column width
			for (size_t i = first_data_idx; i < lines.size(); ++i) {
				std::stringstream ss(lines[i]);
				std::string segment;
				while (std::getline(ss, segment, '\t')) {
					if ((int)segment.length() > col_width)
						col_width =
							(int)segment.length();
				}
			}
			col_width = std::min(col_width + 2, width / 2);
		}

		int label_width = this->pd.lh;
		bool use_frozen_labels =
			(p_id == PaneID::C || p_id == PaneID::D ||
			 p_id == PaneID::E || p_id == PaneID::F);
		if (use_frozen_labels && label_width > safe_w * 0.4) {
			label_width = safe_w * 0.4;
		}

		// --- 3. DATA RENDERING WITH VISUAL OFFSET ---
		// We must find which data index corresponds to our current
		// scroll_offset
		int current_data_idx = first_data_idx;
		int visual_rows_skipped = 0;

		while (current_data_idx < (int)lines.size()) {
			int rows_for_this_item =
				pv_cmp::contains(pd.group_lines,
						 (pt::u32)current_data_idx)
					? 2
					: 1;

			if (visual_rows_skipped + rows_for_this_item >
			    scroll_offset) {
				// This is where we start drawing.
				// If scroll_offset falls exactly on a
				// separator, we handle that offset here.
				break;
			}
			visual_rows_skipped += rows_for_this_item;
			current_data_idx++;
		}

		int y_occupied = 0;
		for (int i = current_data_idx;
		     i < (int)lines.size() && y_occupied < body_h; ++i) {

			// A. Handle Separator
			if (pv_cmp::contains(pd.group_lines, (pt::u32)i)) {
				// If the scroll offset is inside the
				// separator/line pair, we might skip the
				// separator
				if (y_occupied == 0 &&
				    visual_rows_skipped < scroll_offset) {
					// Skip drawing separator because it's
					// scrolled off
				}
				else {
					wattron(win, COLOR_PAIR(3));
					mvwhline(win, data_start_y + y_occupied,
						 1, ACS_HLINE, safe_w);
					wattroff(win, COLOR_PAIR(3));
					y_occupied++;
				}
			}

			if (y_occupied >= body_h)
				break;

			// B. Handle Data Line
			int visual_y = data_start_y + y_occupied;
			bool is_selected = (is_focused && i == selected_line);

			// Clear line and apply selection attribute
			wattrset(win, A_NORMAL);
			if (is_selected)
				wattron(win, A_REVERSE | COLOR_PAIR(1));
			else
				wattrset(win, A_NORMAL);
			mvwprintw(win, visual_y, 1, "%*s", safe_w,
				  ""); // Background fill

			if (use_frozen_labels) {
				// Logic for split panes (Label | Content)
				const line_metadata &meta = pd.meta[i];
				std::string row_label =
					lines[i].substr(0, meta.ref_name_pos);
				std::string row_content =
					lines[i].substr(meta.ref_name_pos);

				// Draw Label
				if (!is_selected)
					wattron(win, A_BOLD | COLOR_PAIR(2));
				else
					wattron(win, A_BOLD);

				mvwprintw(win, visual_y, 1, "%*.*s",
					  label_width, label_width,
					  row_label.c_str());

				if (!is_selected)
					wattroff(win, A_BOLD | COLOR_PAIR(2));
				else
					wattroff(win, A_BOLD);

				// wattroff(win, A_BOLD | COLOR_PAIR(2));

				mvwaddch(win, visual_y, 1 + label_width,
					 ACS_VLINE | COLOR_PAIR(2));

				// Draw Content
				int content_x = 1 + label_width + 2;
				int content_w = safe_w - label_width - 2;
				std::string display_str =
					(row_content.length() >
					 (size_t)horiz_offset)
						? row_content.substr(
							  horiz_offset,
							  content_w)
						: "";

				if (sep_cols) {
					draw_tabular_line(visual_y, content_x,
							  display_str,
							  col_width);
				}
				else {
					for (int j = 0;
					     j < (int)display_str.length();
					     ++j) {
						pt::u32 global_char_pos =
							meta.ref_name_pos +
							horiz_offset + j;
						bool dim = true;
						for (const auto &slice :
						     meta.at_str_slices) {
							if (global_char_pos >=
								    slice.start() &&
							    global_char_pos <
								    (slice.start() +
								     slice.len())) {
								dim = false;
								break;
							}
						}
						// Only dim if we aren't trying
						// to highlight the line
						if (dim & !is_selected)
							wattron(win, A_DIM);
						mvwaddch(win, visual_y,
							 content_x + j,
							 display_str[j]);
						wattroff(win, A_DIM);
					}
				}
			}
			else {
				// Standard logic (Pane A/B)
				std::string display_str =
					(lines[i].length() >
					 (size_t)horiz_offset)
						? lines[i].substr(horiz_offset,
								  safe_w)
						: "";

				bool is_special = pv_cmp::contains(
					pd.special_lines, (pt::u32)i);

				if (is_special && !is_selected)
					wattron(win, COLOR_PAIR(4));

				if (sep_cols)
					draw_tabular_line(visual_y, 1,
							  display_str,
							  col_width);
				else
					mvwprintw(win, visual_y, 1, "%-*.*s",
						  safe_w, safe_w,
						  display_str.c_str());

				wattroff(win, COLOR_PAIR(4));
			}

			if (is_selected)
				wattrset(win, A_NORMAL);

			y_occupied++;
			visual_rows_skipped++; // This helps sync with next loop
					       // iteration
		}

		box(win, 0, 0);
		wrefresh(win);
	}
};

struct pane_params {
	int x;
	int y;
	int width;
	int height;
	PaneID p_id;
	Pane &p;
	bool has_header = false;
	bool sep_cols = false;
	pt::u32 selected_line = 0;
};

// status bar or modeline
// a status bar takes one full line which is saved in the linum field
struct status_bar {
	// --------------------------
	// member variables or fields
	// --------------------------

	WINDOW *win;
	int linum;
	int width;

	// --------------------------
	// constructors & factory fns
	// --------------------------

	status_bar(WINDOW *w, int line_num, int wth)
	    : win(w), linum(line_num), width(wth)
	{}

	// -------
	// methods
	// -------

	void draw(const tui::state::ui_state &state, const std::string &left,
		  const std::string &mid, const std::string &right)
	{
		move(this->linum, 0);
		clrtoeol();

		// 1. Determine the Left text (The Prompt)
		std::string final_left = left;
		if (state.current_mode == Mode::SEARCH) {
			final_left = "/" + state.search_query;
		}
		else if (state.current_mode == Mode::JUMP) {
			// Add this flag to your global state or pane
			final_left = "Jump to line: " + state.jump_query;
		}
		else if (state.current_mode == Mode::COMMAND) {
			final_left = ":" + state.command_prompt;
		}

		// 2. Determine final Right text (Calculate this BEFORE
		// calculating positions)
		std::string final_right = right;
		if (state.current_result_idx != -1 &&
		    !state.search_results.empty()) {
			final_right =
				"[Match " +
				std::to_string(state.current_result_idx + 1) +
				"/" +
				std::to_string(state.search_results.size()) +
				"] " + right;
		}

		// 3. Calculate positions based on final lengths
		int mid_pos = (width - (int)mid.length()) / 2;
		int right_pos = width - (int)final_right.length();

		// 4. Draw Left
		// if (current_mode == Mode::SEARCH || current_mode ==
		// Mode::JUMP)
		mvprintw(this->linum, 0, "%s", final_left.c_str());

		// 5. Draw Middle (only if it fits)
		std::string final_mid = pv_cmp::format(
			"{} [{}]", mid, mode_names[state.current_mode]);

		if (mid_pos > (int)final_left.length() + 1)
			mvprintw(this->linum, mid_pos, "%s", final_mid.c_str());

		// 6. Draw Right (only if it doesn't overlap middle)
		if (right_pos > mid_pos + (int)mid.length() + 1) {
			mvprintw(this->linum, right_pos, "%s",
				 final_right.c_str());
		}

		// 7. Cursor Management
		if (state.current_mode == Mode::SEARCH) {
			curs_set(1);
			move(this->linum, (int)final_left.length());
		}
		else {
			curs_set(0);
			move(0, 0);
		}
	}
};

} // namespace zien::components

#endif // ZIEN_COMPONENTS_HPP
