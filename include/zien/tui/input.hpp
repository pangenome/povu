#include <ncurses.h>

#include <liteseq/refs.h> // for ref_walk, ref

#include "zien/components/components.hpp" // for status_bar
#include "zien/tui/state.hpp"		  // for Mode
#include "zien/tui/tui.hpp"

namespace zien::tui::input
{
using namespace zien::components;
using namespace zien::tui::state;

void perform_search(Pane *pane, ui_state &state);
void clear_search(ui_state &state);
// handle search input
void handle_search_input(int ch, Pane *ap, ui_state &state);
void handle_jump_input(int ch, Pane *active, ui_state &state);
void commands_handler(ui_state &state);
// handle search input
void handle_command_input(int ch, Pane *ap, ui_state &state);
void handle_special_states(int ch, tui_context &tc, ui_state &state);
// navigate and initialize search
void nav(int ch, tui_context &tc, ui_state &state);
void handle_normal_state(int ch, tui_context &tc, ui_state &state);

void handle_input(int ch, tui_context &tc, ui_state &state);
}; // namespace zien::tui::input
