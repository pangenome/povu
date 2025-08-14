#ifndef PV_SUBCOMMANDS_INF_HPP
#define PV_SUBCOMMANDS_INF_HPP

#include "../../include/common/types/compat.hpp"
#include "../../include/graph/bidirected.hpp"
#include "../cli/app.hpp"
#include "./common.hpp"

namespace povu::subcommands::info {
namespace bd = povu::bidirected;
namespace pt = povu::types;

void do_info(const core::config &app_config);
}

#endif
