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
#include "zien/components/paths.hpp"	  // for update_paths
#include "zien/components/refs.hpp"	  // for foo
#include "zien/components/repeats.hpp"	  // for update_repeats
#include "zien/components/vcf.hpp"	  // for comp_vcfs
#include "zien/tui/input.hpp"		  // for
#include "zien/tui/state.hpp"		  // for Mode

namespace zien::tui
{
namespace lq = liteseq;

using namespace zien::tui::state;
using namespace zien::tui::input;
using namespace zien::components;

const pt::u32 REPEATS_PANE_COUNT = 5;
const pt::u32 NO_REPEATS_PANE_COUNT = 4;

const pt::u32 INITIAL_VCF_REC_IDX = 0; // skip header line

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

void create_views(const ui_state &state, Pane &a, Pane &b, Pane &c, Pane &d,
		  Pane *e, Pane &f)
{
	int screen_h = state.screen_h;
	int screen_w = state.screen_w;

	/*
	  --------------
	  Variation view
	  --------------
	*/

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

	/*
	  --------------
	  Paths view
	  --------------
	*/
	setup_pane({0, 0, screen_w, screen_h - 1, PaneID::F, f, false, false});
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

void update_status_bar(ui_state &state, status_bar &sb)
{
	// left
	//
	//
	// middle
	// const std::string &mode = mode_names.at(current_mode);
	const std::string &pane_name = pane_names.at(state.active_pane_id);

	std::string vcf_rec_line; // right
	if (state.current_view == View::PATHS) {
		auto [a, b] = state.paths_view_range;
		vcf_rec_line = pv_cmp::format("[{}/{}]", a, b);
	}
	else {
		// vcf_rec_line =
		//	pv_cmp::format("[{}/{}]", state.vcf_selected_rec + 1,
		//		       vcf_file.record_count());
	}

	sb.draw(state, "", pane_name, vcf_rec_line);
};

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

void draw_all_panes(const bd::VG &g, const mto::from_vcf::VCFile &vcf_file,
		    Pane &top_left_pane, Pane &top_right_pane,
		    Pane &bottom_left_pane, Pane &bottom_right_pane,
		    Pane &repeats_pane, Pane &paths_pane, ui_state &state,
		    status_bar &sb, bool initial_draw)
{

	if (initial_draw) {
		// Start at first record (skip header)
		top_left_pane.selected_line = 1;

		// rest
		top_right_pane.selected_line = 1;
		bottom_left_pane.selected_line = 0;
		bottom_right_pane.selected_line = 0;
	}

	if (state.current_view == View::PATHS) {
		// setup_paths_window(state, paths_pane);
		paths_pane.draw(state);
		update_status_bar(vcf_file, state, sb);
		return;
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

void update_paths_view(const bd::VG &g, ui_state &state, display_lines &pd)
{
	if (state.update_paths_view == false)
		return;

	pd.reset();
	zien::components::paths::update_paths(g, state, pd);
	state.update_paths_view = false;
}

void update_paths_view(const bd::VG &g, ui_state &state, display_lines &pd,
		       status_bar &sb)
{
	if (state.update_paths_view == false)
		return;

	pd.reset();
	zien::components::paths::update_paths(g, state, pd);
	update_status_bar(state, sb);
	state.update_paths_view = false;
}

void update_variants_view(const bd::VG &g,
			  const mto::from_vcf::VCFile &vcf_file,
			  Pane &top_right_pane, Pane &bottom_left_pane,
			  Pane &bottom_right_pane, ui_state &state)
{
	if (state.active_pane_id != PaneID::A)
		return;

	top_right_pane.pd.lines.clear();
	zien::components::genotypes::update_haps(
		g, vcf_file, state.vcf_selected_rec, top_right_pane.pd);

	bottom_left_pane.pd.lines.clear();
	zien::components::refs::update_refs(g, vcf_file, state.vcf_selected_rec,
					    bottom_left_pane.pd);

	bottom_right_pane.pd.lines.clear();
	zien::components::alts::update_alts(g, vcf_file, state.vcf_selected_rec,
					    bottom_right_pane.pd);
}

void update_depenent_panes(const bd::VG &g,
			   const mto::from_vcf::VCFile &vcf_file,
			   Pane &top_right_pane, Pane &bottom_left_pane,
			   Pane &bottom_right_pane, Pane &paths_pane,
			   ui_state &state)
{
	switch (state.current_view) {
	case View::VARIATION:
		update_variants_view(g, vcf_file, top_right_pane,
				     bottom_left_pane, bottom_right_pane,
				     state);
		return;
	case View::PATHS:
		update_paths_view(g, state, paths_pane.pd);
		return;
	}
}

void view_gfa(const bd::VG &g)
{
	tui_context tc;
	ui_state &state = tc.get_state();

	Pane *a = (tc.get_pane(PaneID::A));
	Pane *b = (tc.get_pane(PaneID::B));
	Pane *c = (tc.get_pane(PaneID::C));
	Pane *d = (tc.get_pane(PaneID::D));
	// Pane *e = (tc.get_pane(PaneID::E));
	Pane *f = (tc.get_pane(PaneID::F));

	create_views(state, *a, *b, *c, *d, nullptr, *f);

	status_bar sb(stdscr, state.screen_h - 1, state.screen_w);

	state.current_view = zien::tui::state::View::PATHS;
	state.active_pane_id = PaneID::F; // Initial focus
	state.update_paths_view = true;

	update_paths_view(g, state, f->pd, sb);
	f->draw(state);
	refresh(); // Push all changes to the physical terminal

	int ch{};
	while ((ch = getch()) != 'q') {
		handle_input(ch, tc, state);

		// Updates dependent panes if necessary
		update_paths_view(g, state, f->pd, sb);
		f->draw(state);
		update_status_bar(state, sb);

		refresh(); // Push all changes to the physical terminal
	};

	return;
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
	Pane *f = (tc.get_pane(PaneID::F));

	create_views(state, *a, *b, *c, *d, nullptr, *f);

	status_bar sb(stdscr, state.screen_h - 1, state.screen_w);

	state.active_pane_id = PaneID::A; // Initial focus

	/* top left pane */
	// called once at the start
	zien::components::vcf::comp_vcfs(g, vcf_file, invalid_recs, a->pd);

	// 2. Explicitly trigger the data population BEFORE the loop
	b->pd.lines.clear();
	zien::components::genotypes::update_haps(g, vcf_file,
						 INITIAL_VCF_REC_IDX, b->pd);
	c->pd.lines.clear();
	zien::components::refs::update_refs(g, vcf_file, INITIAL_VCF_REC_IDX,
					    c->pd);
	d->pd.lines.clear();
	zien::components::alts::update_alts(g, vcf_file, INITIAL_VCF_REC_IDX,
					    d->pd);

	d->pd.lines.clear();
	zien::components::paths::update_paths(g, state, f->pd);

	draw_all_panes(g, vcf_file, *a, *b, *c, *d, *e, *f, state, sb, true);

	// aka has repeats
	// does the current record cover a repetitive region?
	bool is_tangled{false};

	int ch{};
	while ((ch = getch()) != 'q') {
		handle_input(ch, tc, state);

		// Updates dependent panes if necessary
		update_depenent_panes(g, vcf_file, *b, *c, *d, *f, state);

		// 2. Tangle Logic (Check if we need to change layout)
		is_tangled = vcf_file.get_records()
				     .at(state.vcf_selected_rec)
				     .is_tangled();

		state.pane_count =
			is_tangled ? REPEATS_PANE_COUNT : NO_REPEATS_PANE_COUNT;

		if (state.toggle_repeats_pane != is_tangled) {
			create_views(state, *a, *b, *c, *d,
				     (is_tangled ? e : nullptr), *f);
			state.toggle_repeats_pane = is_tangled;

			erase();
			refresh();
		}

		draw_all_panes(g, vcf_file, *a, *b, *c, *d, *e, *f, state, sb,
			       false);
		refresh(); // Push all changes to the physical terminal
	};

	return;
}

} // namespace zien::tui
