// ----------------------------------------------------------------------------
// Simple ncurses-based dual-pane viewer with horizontal and vertical scrolling
// ----------------------------------------------------------------------------

#include "zien/tui/tui.hpp"

// std
#include <sstream>
#include <string>
#include <vector>

#include <atomic>
// #include <thread>
#include <vector>

// ncurses
#include <ncurses.h>

// liteseq
#include <liteseq/refs.h> // for ref_walk, ref

// povu
#include "mto/from_vcf.hpp" // for VCFile

#include "povu/common/core.hpp"	     // for pt
#include "povu/common/utils.hpp"     // for pu::concat_with
#include "povu/graph/bidirected.hpp" // for VG

namespace zien::tui
{
namespace lq = liteseq;

// Use Group Separator (0x1D) as the invisible delimiter
// const char SEP = '\x1D';

std::string search_query = "";
bool search_mode = false;
std::vector<int> search_results;
int current_result_idx = -1;

enum class PaneID : uint8_t {
	A, // top left
	B, // top right
	C, // bottom left
	D  // bottom right
};

std::map<PaneID, std::string> pane_names = {{PaneID::A, "VCF"},
					    {PaneID::B, "Haps"},
					    {PaneID::C, "Refs"},
					    {PaneID::D, "Alts"}};

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
		// Height - 2 (Borders). If header exists, subtract 2 more (Text
		// + HLine).
		int data_area_h = height - 2 - (has_header ? 2 : 0);
		int first_data_idx = has_header ? 1 : 0;

		// Relative position of selection compared to the first data
		// line
		int relative_pos = selected_line - first_data_idx;

		if (relative_pos < scroll_offset)
			scroll_offset = relative_pos;

		if (relative_pos >= scroll_offset + data_area_h)
			scroll_offset = relative_pos - data_area_h + 1;

		if (scroll_offset < 0)
			scroll_offset = 0;
	}

	void draw_tabular_line(WINDOW *win, int y, int x,
			       const std::string &line, int col_width)
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
		pt::u32 max_idx = (pt::u32)pd.lines.size() - 1;
		if ((pt::u32)selected_line < max_idx) {
			selected_line++;
		}
	}

	void scroll_up()
	{
		pt::u32 min_idx = has_header ? 1 : 0;
		if ((pt::u32)selected_line > min_idx) {
			selected_line--;
		}
	}

	void draw(bool is_focused)
	{
		// int col_width = 30;
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
			(p_id == PaneID::C || p_id == PaneID::D);

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

			if (use_frozen_labels &&
			    pv_cmp::contains(this->pd.group_lines,
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
					draw_tabular_line(
						win, visual_y, content_x,
						display_str, col_width);
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
					draw_tabular_line(win, visual_y, 1,
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

void perform_search(Pane *pane, const std::string &query)
{
	search_results.clear();
	for (pt::u32 i = 0; i < pane->pd.lines.size(); ++i)
		if (pane->pd.lines[i].find(query) != std::string::npos)
			search_results.push_back(i);

	if (!search_results.empty()) {
		current_result_idx = 0;
		pane->selected_line = search_results[0];
	}
}

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
	p.selected_line = sl;
}

void setup_layout(int screen_h, int screen_w, Pane &a, Pane &b, Pane &c,
		  Pane &d)
{
	int mid_h = screen_h / 2;
	int mid_w = screen_w / 2;

	int gutter = 1; // gutter of 1 char between windows

	setup_pane({0, 0, mid_w - gutter, mid_h, PaneID::A, a, true, false, 1});
	setup_pane({mid_w, 0, screen_w - mid_w, mid_h, PaneID::B, b, false,
		    true, 1});

	// 2. Bottom Left: Half width, bottom half
	int bottom_h = screen_h - mid_h;
	// We subtract 1 from bottom_h if we want to leave the very last
	// line for the search bar
	int actual_view_h = bottom_h - 1;

	setup_pane({0, mid_h, mid_w - gutter, actual_view_h, PaneID::C, c});

	// 3. Bottom Right: Remaining width, bottom half
	setup_pane(
		{mid_w, mid_h, screen_w - mid_w, actual_view_h, PaneID::D, d});
}

char lq_strand_to_or_e(lq::strand s)
{
	return (s == lq::strand::STRAND_FWD) ? '>' : '<';
}

constexpr povu::refs::ref_format_e PN = povu::refs::ref_format_e::PANSN;

std::set<pt::id_t> get_ref_ids(const bd::VG &g, const std::string &sn,
			       pt::u32 phase_idx)
{
	std::set<pt::id_t> filtered_ref_ids;

	for (pt::id_t r_id : g.get_refs_in_sample(sn)) {
		const povu::refs::Ref &r = g.get_ref_by_id(r_id);

		// TODO: find a better way to handle non PANSN
		if (r.get_format() != PN) // just trust it
			filtered_ref_ids.insert(r_id);
		else if (r.get_format() == PN && r.get_hap_id() == phase_idx)
			filtered_ref_ids.insert(r_id);
	}

	return filtered_ref_ids;
}

pt::id_t extract_anchor_v_id(const std::string &at)
{
	std::string v_id_str = "";
	// zero is a > or <
	for (pt::u32 i{1}; i < at.size(); i++) {
		char c = at[i];
		if (c == '>' || c == '<')
			break;

		if (std::isdigit(c))
			v_id_str += c;
	}

	return static_cast<pt::id_t>(std::stoll(v_id_str));
}

pt::u32 at_str_step_count(const std::string &at_str)
{
	pt::u32 count{};
	for (pt::u32 i{}; i < at_str.size(); i++) {
		char c = at_str[i];
		if (c == '>' || c == '<')
			count++;
	}

	return count;
}

void comp_update_refs(const bd::VG &g, const mto::from_vcf::VCFile &vcf_file,
		      pt::u32 selected_rec_idx, pt::u32 at_idx,
		      display_lines &pd)
{
	pt::u32 line_count = pd.lines.size();
	pt::u32 ctr{line_count}; // counter for line metadata

	// what we compute and return
	// lm &meta = pd.meta;
	std::vector<std::string> &ref_lines = pd.lines;

	const mto::from_vcf::VCFRecord &rec =
		vcf_file.get_records().at(selected_rec_idx);
	const mto::from_vcf::gt_data &d = rec.get_genotypes();
	const std::vector<mto::from_vcf::at_meta> at_meta =
		d.get_data().at(at_idx);

	const std::string &s = rec.get_at(at_idx);
	pt::u32 at_str_len = s.length();
	pt::u32 at_str_sc = at_str_step_count(s);
	pt::id_t anchor_v_id = extract_anchor_v_id(s);

	for (const auto &[sample_idx, phase_idx] : at_meta) {
		std::string sn = vcf_file.get_sample_name(sample_idx);
		std::string curr_l = "";

		for (pt::u32 h_idx : get_ref_ids(g, sn, phase_idx + 1)) {

			const lq::ref_walk *rw = g.get_ref_vec(h_idx)->walk;
			pt::u32 N = rw->step_count;

			curr_l.clear();

			// curr_l += "[" + std::to_string(N) + "] ";
			curr_l += g.get_tag(h_idx);

			pt::u32 header_len = curr_l.length();

			if (header_len > pd.lh)
				pd.lh = header_len;

			// std::cerr << "Anchor v_id: " << anchor_v_id << "\t"
			//	  << sn << "\n";

			pt::u32 anchor_v_idx = g.v_id_to_idx(anchor_v_id);
			const std::vector<pt::idx_t> &starts =
				g.get_vertex_ref_idxs(anchor_v_idx, h_idx);

			if (starts.empty())
				continue;

			line_metadata &lm = pd.meta[ctr]; // create if
							  // not exists
			lm.ref_name_pos = header_len;

			for (pt::u32 s{}; s < starts.size(); s++) {
				pt::u32 ref_pos = starts[s];
				pt::u32 i = ref_pos > 10 ? ref_pos - 10 : 0;
				pt::u32 end =
					std::min(ref_pos + at_str_sc + 10, N);

				for (; i < end; i++) {

					if (i == ref_pos) {
						pt::slice w{
							(pt::u32)curr_l.size(),
							at_str_len};

						lm.at_str_slices.emplace_back(
							w);
					}

					pt::id_t v_id = rw->v_ids[i];
					char o = lq_strand_to_or_e(
						rw->strands[i]);

					curr_l += o;
					curr_l += std::to_string(v_id);
				}

				if (s + 1 < starts.size())
					curr_l += " ... ";

				// curr_l += "\t";
			}

			ref_lines.push_back(curr_l);
			ctr++;
			// pd.lines.push_back(curr_l);
			//  pd.p->lines.push_back(curr_l);
		}
	}

	return;
};

std::vector<std::string> comp_vcf_lines(const mto::from_vcf::VCFile &vcf_file)
{
	std::vector<std::string> lines;
	std::stringstream ss; // string buffer

	// header
	const std::vector<std::string> COL_NAMES = {
		"VAR TYPE", "POS", "ID", "REF", "ALT", "REF AT", "ALT AT"};

	lines.emplace_back(pu::concat_with(COL_NAMES, '\t'));

	const std::vector<mto::from_vcf::VCFRecord> &all_recs =
		vcf_file.get_records();
	for (const mto::from_vcf::VCFRecord &rec : all_recs) {
		rec.viewer(ss);
		lines.push_back(ss.str());
		ss.str(std::string()); // Clear the stringstream
	}

	return lines;
}

std::vector<std::string> comp_gt_data(const mto::from_vcf::VCFile &vcf_file,
				      pt::u32 selected_rec_idx)
{
	const mto::from_vcf::VCFRecord &rec =
		vcf_file.get_records().at(selected_rec_idx);
	const mto::from_vcf::gt_data &d = rec.get_genotypes();

	pt::u32 at_count = rec.get_at_count();

	std::vector<std::string> hap_lines;
	std::string hl;

	for (pt::u32 at_idx{}; at_idx < at_count; at_idx++) {
		const std::string &at = rec.get_at(at_idx);
		hl += at;
		hl += "\t";
	}

	hap_lines.emplace_back(hl);
	hl.clear();

	// count rows
	pt::u32 row_count{};
	for (pt::u32 at_idx{}; at_idx < at_count; at_idx++) {
		if (!pv_cmp::contains(d.get_data(), at_idx)) {
			continue;
		}
		// try {
		//	const std::vector<povu::io::from_vcf::at_meta> at_meta =
		//		d.get_data().at(at_idx);
		// }
		// catch (const std::out_of_range &e) {
		//	continue;
		// }

		const std::vector<mto::from_vcf::at_meta> &at_meta =
			d.get_data().at(at_idx);
		pt::u32 N = at_meta.size();
		if (N > row_count)
			row_count = N;
	}

	for (pt::u32 row_idx{}; row_idx < row_count; row_idx++) {
		for (pt::u32 at_idx{}; at_idx < at_count; at_idx++) {

			if (!pv_cmp::contains(d.get_data(), at_idx)) {
				// AT is not supported by any haps
				hl += "?";
				hl += "\t";
				continue;
			}

			// try {
			//	const std::vector<povu::io::from_vcf::at_meta>
			//		at_meta = d.get_data().at(at_idx);
			// }
			// catch (const std::out_of_range &e) {
			//	hl += "\t";
			//	continue;
			// }

			const std::vector<mto::from_vcf::at_meta> &at_meta =
				d.get_data().at(at_idx);

			if (row_idx < at_meta.size()) {
				const auto &[sample_idx, phase_idx] =
					at_meta[row_idx];
				std::string sn =
					vcf_file.get_sample_name(sample_idx);
				std::string l = sn + '#' +
						std::to_string(phase_idx + 1);
				hl += l;
				hl += "\t";
			}
			else {
				hl += "\t";
			}
		}
		hap_lines.emplace_back(hl);
		hl.clear();
	}

	return hap_lines;
};

// navigate and initialize search
void nav(int ch, Pane *&active, const std::map<PaneID, Pane *> &panes)
{
	PaneID active_p_id = active->p_id;

	switch (ch) {
	case '/':
		search_mode = true;
		search_query = "";
		break;
	case 'n':
		if (!search_results.empty()) {
			current_result_idx = (current_result_idx + 1) %
					     search_results.size();
			active->selected_line =
				search_results[current_result_idx];
		}
		break;
	case '\t': // Cycle: Top -> Left -> Right -> Top
		switch (active_p_id) {
		case PaneID::A:
			active = panes.at(PaneID::B);
			break;
		case PaneID::B:
			active = panes.at(PaneID::C);
			break;
		case PaneID::C:
			active = panes.at(PaneID::D);
			break;
		default: // also covers PaneID::D
			active = panes.at(PaneID::A);
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
}

// handle search input
void handle_search_input(int ch, Pane *&active)
{
	switch (ch) {
	case '\n':
		search_mode = false;
		perform_search(active, search_query);
		break;
	case 27: // ESC key
		search_mode = false;
		break;
	case KEY_BACKSPACE:
	case 127: // Handle Backspace or DEL
		if (!search_query.empty())
			search_query.pop_back();
		break;
	default:
		if (ch >= 32 && ch <= 126) // Only add printable characters
			search_query += (char)ch;
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

		// Fetch dimensions inside the loop in case the user resizes the
		// terminal
		int h, w;
		getmaxyx(stdscr, h, w);

		// No attron() or COLOR_PAIR needed for default white on black
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
	erase();   // Clears the internal buffer (removes the box and text)
	refresh(); // Pushes that empty buffer to the physical screen

	curs_set(1); // restore the cursor

	// Reset ncurses to blocking mode for the rest of the app
	timeout(-1);
}

void view(const bd::VG &g, const mto::from_vcf::VCFile &vcf_file,
	  const std::vector<pt::u32> &invalid_recs)
{
	// ------------------------
	// colours
	//
	// ansi colours lookup
	// https://gist.github.com/JBlond/2fea43a3049b38287e5e9cefc87b2124
	// ------------------------

	// Use indices in the 240+ range to stay safe
	short dark_gray_idx = 234;   // Very dark (near black)
	short medium_gray_idx = 240; // Classic gray

	init_pair(1, COLOR_WHITE, COLOR_BLACK);	    // pure white (Color 7)
	init_pair(2, medium_gray_idx, COLOR_BLACK); // medium gray
	init_pair(3, dark_gray_idx, COLOR_BLACK); // dark gray for the separator
	init_pair(4, COLOR_RED, COLOR_BLACK);	  // bright red for errors

	int screen_h, screen_w;
	getmaxyx(stdscr, screen_h, screen_w);

	Pane top_left_pane, top_right_pane, bottom_left_pane, bottom_right_pane;
	const std::map<PaneID, Pane *> PANES = {
		{PaneID::A, &top_left_pane},
		{PaneID::B, &top_right_pane},
		{PaneID::C, &bottom_left_pane},
		{PaneID::D, &bottom_right_pane}};
	setup_layout(screen_h, screen_w, top_left_pane, top_right_pane,
		     bottom_left_pane, bottom_right_pane);

	// Initial focus
	Pane *active = &top_left_pane;

	/* top left pane */
	std::vector<std::string> vcf_lines = comp_vcf_lines(vcf_file);
	top_left_pane.pd.lines.swap(vcf_lines);
	top_left_pane.pd.special_lines.insert(invalid_recs.begin(),
					      invalid_recs.end());

	/* top-right (haplotype) pane */
	auto update_haps = [&](pt::u32 selected_rec_idx)
	{
		std::vector<std::string> hap_lines =
			comp_gt_data(vcf_file, selected_rec_idx);

		for (const std::string &hl : hap_lines)
			top_right_pane.pd.lines.push_back(hl);
	};

	/* left (ref) pane */
	auto update_refs = [&](pt::u32 selected_rec_idx, display_lines &pd)
	{
		const pt::u32 REF_AT_IDX{0};
		pd.reset();
		comp_update_refs(g, vcf_file, selected_rec_idx, REF_AT_IDX, pd);
	};

	/* right (alt) pane */
	auto update_alts = [&](pt::u32 selected_rec_idx, display_lines &pd)
	{
		pd.reset();

		const mto::from_vcf::VCFRecord &rec =
			vcf_file.get_records().at(selected_rec_idx);

		pt::u32 N = rec.get_alt_count() + 1;

		for (pt::u32 i{1}; i < N; i++) { // start at 1 to skip REF
			comp_update_refs(g, vcf_file, selected_rec_idx, i, pd);

			if (i + 1 < N) // add a separator line between alleles
				pd.group_lines.insert(pd.lines.size());
		}
	};

	auto draw_all_panes = [&](bool initial_draw = false)
	{
		top_left_pane.draw(active == &top_left_pane);
		top_right_pane.draw(active == &top_right_pane);
		bottom_left_pane.draw(active == &bottom_left_pane);
		bottom_right_pane.draw(active == &bottom_right_pane);

		if (initial_draw)
			refresh();
	};

	// 2. Explicitly trigger the data population BEFORE the loop
	top_right_pane.pd.lines.clear();
	update_haps(0);
	bottom_left_pane.pd.lines.clear();
	update_refs(0, bottom_left_pane.pd);
	bottom_right_pane.pd.lines.clear();
	update_alts(0, bottom_right_pane.pd);

	draw_all_panes(true);

	int ch = 0;
	do {
		// 1. Handle Input
		if (search_mode)
			handle_search_input(ch, active);
		else
			nav(ch, active, PANES); // nav and search

		// 2. Update dependent panes if necessary
		if (active == &top_left_pane) {
			pt::u32 vcf_data_row = active->selected_line - 1;

			top_right_pane.pd.lines.clear();
			update_haps(vcf_data_row);

			bottom_left_pane.pd.lines.clear();
			update_refs(vcf_data_row, bottom_left_pane.pd);

			bottom_right_pane.pd.lines.clear();
			update_alts(vcf_data_row, bottom_right_pane.pd);
		}

		draw_all_panes();

		// 3. Draw Search Bar (Always do this)
		move(screen_h - 1, 0); // Move to the very bottom line
		clrtoeol(); // Clear the entire line from cursor to end
		if (search_mode) {
			attron(COLOR_PAIR(1));
			printw("/%s", search_query.c_str());
			attroff(COLOR_PAIR(1));
		}
		else if (!search_query.empty()) {
			printw("Search: %s (%d/%zu)", search_query.c_str(),
			       current_result_idx + 1, search_results.size());
		}

		refresh(); // Push all changes to the physical terminal
	} while ((ch = getch()) != 'q');

	return;
}

} // namespace zien::tui
