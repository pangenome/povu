#ifndef PV_SUBCOMMANDS_DEC_HPP
#define PV_SUBCOMMANDS_DEC_HPP

#include <string_view> // for string_view

#include "povu/common/app.hpp"

namespace povu::subcommands::decompose
{
constexpr std::string_view MODULE = "povu::subcommands::decompose";

void do_decompose(const core::config &app_config);
} // namespace povu::subcommands::decompose

#endif
