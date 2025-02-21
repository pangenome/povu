#ifndef PV_SUBCOMMANDS_INF_HPP
#define PV_SUBCOMMANDS_INF_HPP

#include "../../include/graph/bidirected.hpp"
#include "../cli/app.hpp"
//#include "../io/io.hpp"
#include "./common.hpp"

namespace povu::subcommands::info {
  //using namespace povu::io;
namespace bd = povu::bidirected;

void do_info(const core::config &app_config);
}

#endif
