#ifndef PV_CLI_HPP
#define PV_CLI_HPP

#include <string> // for basic_string, string

#include "povu/common/app.hpp" // for config

namespace cli
{

// app version string constant
const std::string VERSION = "0.0.0-alpha";

int cli(int argc, char **argv, core::config &app_config);
} // namespace cli

#endif // PV_CLI_HPP
