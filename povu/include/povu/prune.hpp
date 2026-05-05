#ifndef PV_SUBCOMMANDS_PRUNE_HPP
#define PV_SUBCOMMANDS_PRUNE_HPP

#include <string_view> // for string_view

#include <quilt/app.hpp> // for config

namespace povu::subcommands::prune
{

void do_prune(const core::config &app_config);
} // namespace povu::subcommands::prune

#endif
