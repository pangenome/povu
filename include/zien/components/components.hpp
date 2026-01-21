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

		// Calculate how many extra visual lines (separators) are above
		// our selection
		int separators_above = 0;
		for (auto line_idx : pd.group_lines) {
			if ((int)line_idx <= selected_line) {
				separators_above++;
			}
		}

		// The virtual Y position is the data index + separators
		int virtual_selected_pos =
			(selected_line - first_data_idx) + separators_above;

		if (virtual_selected_pos < scroll_offset)
			scroll_offset = virtual_selected_pos;

		if (virtual_selected_pos >= scroll_offset + data_area_h)
			scroll_offset = virtual_selected_pos - data_area_h + 1;

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

		update_scroll();
		werase(win);
		box(win, 0, 0);

		const std::vector<std::string> &lines = this->pd.lines;

		if (!title.empty())
			mvwprintw(win, 0, 2, " %s ", title.c_str());

		int safe_w = width - 2;
		int data_start_y = has_header ? 3 : 1;
		int body_h = (height - 2) - (has_header ? 2 : 0);
		int first_data_idx = has_header ? 1 : 0;

		// --- DYNAMIC COL_WIDTH CALCULATION ---
		int dynamic_col_width = 10; // Minimum default
		if (sep_cols) {
			for (int i = 0; i < body_h; i++) {
				int line_idx =
					first_data_idx + scroll_offset + i;
				if (line_idx >= (int)lines.size())
					break;

				std::stringstream ss(lines[line_idx]);
				std::string segment;
				while (std::getline(ss, segment, '\t')) {
					if ((int)segment.length() >
					    dynamic_col_width) {
						dynamic_col_width =
							(int)segment.length();
					}
				}
			}
			dynamic_col_width += 2; // Add gutter for readability

			// Cap col_width to 50% of window so we don't push
			// secondary columns entirely off-screen
			if (dynamic_col_width > (width / 2))
				dynamic_col_width = width / 2;
		}

		// Now use dynamic_col_width instead of the hardcoded 30
		int col_width = dynamic_col_width;

		// --- DYNAMIC LABEL WIDTH CALCULATION ---
		int label_width = this->pd.lh;
		bool use_frozen_labels =
			(p_id == PaneID::C || p_id == PaneID::D ||
			 p_id == PaneID::E);

		if (use_frozen_labels) {
			// for (const auto &line : lines) {
			//	size_t split_pos = line.find_first_of(SEP);
			//	if (split_pos != std::string::npos) {
			//		// The split position is exactly the
			//		// length of the name
			//		if ((int)split_pos > label_width) {
			//			label_width = (int)split_pos;
			//		}
			//	}
			// }
			// Cap label width at 40% of window to preserve space
			// for the path
			// if (label_width > safe_w * 0.4)
			//	label_width = safe_w * 0.4;

			if (this->pd.lh > safe_w * 0.4)
				label_width = safe_w * 0.4;
		}

		// --- HEADER (Pane A/B) ---
		if (has_header && !lines.empty()) {
			wattron(win, A_BOLD | COLOR_PAIR(2));
			mvwprintw(win, 1, 1, "%-*.*s", safe_w, safe_w,
				  lines[0].c_str());
			wattroff(win, A_BOLD | COLOR_PAIR(2));
			mvwhline(win, 2, 1, ACS_HLINE, safe_w);
		}

		// --- DATA ROWS ---
		int separator_count = 0;
		for (int i = 0; i < body_h; i++) {
			int line_idx = first_data_idx + scroll_offset + i;
			if (line_idx >= (int)lines.size())
				break;

			// Shift the actual text row down based on how many
			// separators we've hit
			int visual_y = data_start_y + i + separator_count;

			// Check if we exceed pane height
			if (visual_y >= height - 1)
				break;

			// Draw separator line if needed
			if (pv_cmp::contains(this->pd.group_lines,
					     (pt::u32)line_idx)) {
				wattron(win, COLOR_PAIR(3));
				mvwhline(win, visual_y, 1, ACS_HLINE, safe_w);
				wattroff(win, COLOR_PAIR(3));

				separator_count++;
				visual_y++; // Move the text for THIS line to
					    // the next row

				if (visual_y >= height - 1)
					break;
			}

			// Reset and clear the row background
			wattrset(win, A_NORMAL);
			mvwprintw(win, visual_y, 1, "%*s", safe_w, "");

			bool is_selected =
				(is_focused && line_idx == selected_line);
			if (is_selected)
				wattron(win, A_REVERSE | COLOR_PAIR(1));

			if (use_frozen_labels) {
				int gutter_width = 2;

				const std::string &full_line = lines[line_idx];

				const line_metadata &line_meta =
					this->pd.meta[line_idx];

				pt::u32 split_pos = line_meta.ref_name_pos;

				std::string row_label =
					full_line.substr(0, split_pos);

				std::cerr << __func__ << " Line " << line_idx
					  << " split at " << split_pos << " rl "
					  << row_label << "\n";

				// Draw Label: Using A_BOLD and right-alignment
				// for a cleaner look
				wattron(win, A_BOLD | COLOR_PAIR(2));
				mvwprintw(win, visual_y, 1, "%*.*s",
					  label_width, label_width,
					  row_label.c_str());
				wattroff(win, A_BOLD | COLOR_PAIR(2));

				// Draw Separator (The visual replacement for
				// the colon)
				mvwaddch(win, visual_y, 1 + label_width,
					 ACS_VLINE | COLOR_PAIR(2));

				int content_x = 1 + label_width + gutter_width;
				int content_w =
					safe_w - label_width - gutter_width;

				std::string row_content =
					full_line.substr(split_pos);

				std::string display_str = "";
				if ((int)row_content.length() > horiz_offset)
					display_str = row_content.substr(
						horiz_offset, content_w);

				if (sep_cols) {
					draw_tabular_line(visual_y, content_x,
							  display_str,
							  col_width);
				}
				else {
					for (int j = 0;
					     j < (int)display_str.length();
					     j++) {

						pt::u32 global_idx =
							split_pos +
							horiz_offset + j;

						// 2. Check if global_idx is
						// inside any of the slices
						bool should_dim = true;
						for (const pt::slice &slice :
						     line_meta.at_str_slices) {
							// auto [start, len] =
							//	slice;
							// pt::slice usually has
							// .offset and .len (or
							// .start and .end)
							if (global_idx >=
								    slice.start() &&
							    global_idx <
								    (slice.start() +
								     slice.len())) {
								should_dim =
									false;
								break;
							}
						}

						// 3. Apply attributes
						if (should_dim) {
							wattron(win, A_DIM);
						}
						else {
							wattroff(win, A_DIM);
						}

						mvwaddch(win, visual_y,
							 content_x + j,
							 display_str[j]);
					}

					wattroff(win, A_DIM);
				}
			}
			else {
				// Standard Logic for A & B
				std::string display_str = "";
				if ((int)lines[line_idx].length() >
				    horiz_offset)
					display_str = lines[line_idx].substr(
						horiz_offset, safe_w);

				if (sep_cols) {
					draw_tabular_line(visual_y, 1,
							  display_str,
							  col_width);
				}
				else {
					if (pv_cmp::contains(
						    this->pd.special_lines,
						    (pt::u32)line_idx)) {
						wattron(win, COLOR_PAIR(4));
					}

					mvwprintw(win, visual_y, 1, "%-*.*s",
						  safe_w, safe_w,
						  display_str.c_str());

					wattroff(win, COLOR_PAIR(4));
				}
			}

			// Pad the highlight to fill the remainder of the window
			// width
			int cur_x = getcurx(win);
			if (cur_x < safe_w + 1) {
				whline(win, ' ', (safe_w + 1) - cur_x);
			}

			if (is_selected)
				wattrset(win, A_NORMAL);
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

		// 1. Determine final Left text
		// std::string final_left =
		//	search_mode ? ("/" + search_query) : left;

		// 1. Determine the Left text (The Prompt)
		std::string final_left = left;
		if (state.current_mode == Mode::SEARCH) {
			final_left = "/" + state.search_query;
		}
		else if (state.current_mode == Mode::JUMP) {
			// Add this flag to your global state or pane
			final_left = "Jump to line: " + state.jump_query;
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
