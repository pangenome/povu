#ifndef PV_SUBCOMMANDS_DEC_HPP
#define PV_SUBCOMMANDS_DEC_HPP

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <sys/types.h>
#include <thread>
#include <utility>

#include "../../include/algorithms/concealed.hpp"
#include "../../include/algorithms/flubbles.hpp"
#include "../../include/algorithms/midi.hpp"
#include "../../include/algorithms/parallel.hpp"
#include "../../include/algorithms/smothered.hpp"
#include "../../include/algorithms/tiny.hpp"
#include "../../include/common/tree_utils.hpp"
#include "../../include/common/types/compat.hpp"
#include "../../include/graph/bidirected.hpp"
#include "../../include/graph/spanning_tree.hpp"
#include "../../include/io/pvst.hpp"
#include "../../include/io/to_vcf.hpp"
#include "../cli/app.hpp"
#include "../cli/cli.hpp"
#include "./common.hpp"

namespace povu::subcommands::decompose {
constexpr std::string_view MODULE = "povu::subcommands::decompose";

namespace fs = std::filesystem;

namespace pvt = povu::types::genomics;
namespace pt = povu::types;
namespace bd = povu::bidirected;
namespace pgt = povu::types::graph;
namespace pst = povu::spanning_tree;
namespace pvtr = povu::tree;
namespace pvst = povu::types::pvst;
namespace pic = povu::io::common;
namespace ptu = povu::tree_utils;
namespace pfl = povu::flubbles;

using namespace povu::subcommands::common;

void do_decompose(const core::config &app_config);
}

#endif
