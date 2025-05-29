#ifndef PV_SUBCOMMANDS_DEC_HPP
#define PV_SUBCOMMANDS_DEC_HPP

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <sys/types.h>
#include <thread>
#include <utility>

#include "../../include/algorithms/hubbles.hpp"
#include "../../include/algorithms/algorithms.hpp"
#include "../../include/common/types.hpp"
#include "../../include/graph/pvst.hpp"
#include "../../include/graph/spanning_tree.hpp"
#include "../../include/graph/bidirected.hpp"
#include "../io/to_vcf.hpp"
#include "../io/pvst.hpp"
#include "../cli/app.hpp"
#include "../cli/cli.hpp"
#include "./common.hpp"

namespace povu::subcommands::deconstruct {
namespace fs = std::filesystem;

namespace pvt = povu::types::genomics;
namespace pt = povu::types;
namespace bd = povu::bidirected;
namespace pgt = povu::graph_types;
namespace pst = povu::spanning_tree;
namespace pgt = povu::graph_types;
namespace pvtr = povu::tree;
namespace pvst = povu::types::pvst;
namespace pic = povu::io::common;


using namespace povu::subcommands::common;

void do_deconstruct(const core::config &app_config);
}

#endif
