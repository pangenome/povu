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
#include "../../include/common/compat.hpp"
#include "../../include/graph/bidirected.hpp"
#include "../../include/graph/spanning_tree.hpp"
#include "../../include/graph/tree_utils.hpp"
#include "../../include/io/from_gfa.hpp"
#include "../../include/io/to_pvst.hpp"
#include "../../include/io/to_vcf.hpp"
#include "common/app.hpp"
#include "../cli/cli.hpp"

namespace povu::subcommands::decompose
{
constexpr std::string_view MODULE = "povu::subcommands::decompose";

namespace fs = std::filesystem;
namespace bd = povu::bidirected;
namespace pgt = povu::types::graph;
namespace pst = povu::spanning_tree;
namespace pvst = povu::pvst;
namespace pic = povu::io::common;
namespace ptu = povu::tree_utils;
namespace pfl = povu::flubbles;

void do_decompose(const core::config &app_config);
} // namespace povu::subcommands::decompose

#endif
