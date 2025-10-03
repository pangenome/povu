#ifndef PV_SUBCOMMANDS_INF_HPP
#define PV_SUBCOMMANDS_INF_HPP

#include "../../include/common/compat.hpp"
#include "../../include/graph/bidirected.hpp"
#include "../../include/io/from_gfa.hpp"
#include "common/app.hpp"

namespace povu::subcommands::info
{
namespace bd = povu::bidirected;

void do_info(const core::config &app_config);
} // namespace povu::subcommands::info

#endif
