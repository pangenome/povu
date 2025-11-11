#ifndef PV_CLI_HPP
#define PV_CLI_HPP

#include <string_view> // for string_view

#include "povu/common/app.hpp" // for config

namespace cli
{
constexpr std::string_view VERSION = "0.0.1-alpha"; // app & lib version

int cli(int argc, char **argv, core::config &app_config);
} // namespace cli

#endif // PV_CLI_HPP
