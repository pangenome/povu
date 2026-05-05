#ifndef PV_SUBCOMMANDS_VIEW_HPP
#define PV_SUBCOMMANDS_VIEW_HPP

#include <string_view> // for string_view

#include <quilt/app.hpp> // for config

namespace povu::subcommands::view
{

void do_view(const core::config &app_config);
} // namespace povu::subcommands::view

#endif
