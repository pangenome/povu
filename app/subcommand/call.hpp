#ifndef PV_SUBCOMMANDS_CALL_HPP
#define PV_SUBCOMMANDS_CALL_HPP

#include <string_view> // for string_view

#include "povu/common/app.hpp" // for config

namespace povu::subcommands::call
{
constexpr std::string_view MODULE = "povu::subcommands::call";

void do_call(core::config &app_config);
} // namespace povu::subcommands::call
#endif
