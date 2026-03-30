#ifndef PV_SUBCOMMANDS_PRUNE_HPP
#define PV_SUBCOMMANDS_PRUNE_HPP

#include <string_view> // for string_view

#include "povu/common/app.hpp"

namespace povu::subcommands::prune
{
constexpr std::string_view MODULE = "povu::subcommands::prune";

void do_prune(const core::config &app_config);
} // namespace povu::subcommands::prune

#endif
