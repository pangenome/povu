#ifndef PV_PROGRESS_HPP
#define PV_PROGRESS_HPP

#include <cstddef>
#include <indicators/dynamic_progress.hpp>
#include <indicators/indeterminate_progress_bar.hpp>
#include <indicators/multi_progress.hpp>
#include <indicators/progress_bar.hpp>
#include <indicators/setting.hpp>
#include <optional>

namespace povu::progress
{
using namespace indicators;

inline void
set_progress_bar_common_opts(ProgressBar *bar,
			     std::optional<std::size_t> opt_max = std::nullopt)
{
	bar->set_option(option::BarWidth{50});
	bar->set_option(option::Start{"["});
	bar->set_option(option::Fill{"="});
	bar->set_option(option::Lead{">"});
	bar->set_option(option::Remainder{" "});
	bar->set_option(option::End{"]"});
	bar->set_option(option::ForegroundColor{Color::green});
	bar->set_option(option::ShowPercentage{true});
	bar->set_option(option::ShowElapsedTime{true});

	if (opt_max)
		bar->set_option(option::MaxProgress{*opt_max});
}

// progress bar for indeterminate tasks
inline void set_progress_bar_ind(IndeterminateProgressBar *bar)
{
	bar->set_option(option::BarWidth{50});
	bar->set_option(option::Start{"["});
	bar->set_option(option::Fill{"·"});
	bar->set_option(option::Lead{"<==>"});
	bar->set_option(option::End{"]"});
	bar->set_option(option::ForegroundColor{Color::green});
}
}; // namespace povu::progress

#endif
