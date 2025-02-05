#ifndef POVU_SUBCOMMANDS_HPP
#define POVU_SUBCOMMANDS_HPP

#include <thread> // for deconstruct

#include "../../include/algorithms/algorithms.hpp"
#include "../../include/common/types.hpp"
#include "../../include/graph/biedged.hpp"
#include "../../include/graph/flubble_tree.hpp"
#include "../../include/graph/spanning_tree.hpp"
#include "../io/io.hpp"
#include "../cli/app.hpp"
#include "../cli/cli.hpp"


namespace povu::subcommands {
namespace pt = povu::types;
namespace bd = povu::bidirected;
namespace pgt = povu::graph_types;
namespace pst = povu::spanning_tree;
namespace pgt = povu::graph_types;
namespace pvtr = povu::tree;

void do_call(const core::config &app_config);
void do_info(const core::config &app_config);
void do_deconstruct(const core::config &app_config);
}

#endif
