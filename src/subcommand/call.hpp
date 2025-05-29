#ifndef PV_SUBCOMMANDS_CALL_HPP
#define PV_SUBCOMMANDS_CALL_HPP

#include <thread>

#include "../../include/graph/bidirected.hpp"
#include "../io/to_vcf.hpp"

#include "../../include/common/genomics.hpp"
#include "../../include/genomics/genomics.hpp"
#include "../cli/app.hpp"
#include "../cli/cli.hpp"
#include "../io/pvst.hpp"
#include "./common.hpp"

namespace povu::subcommands::call {
namespace fs = std::filesystem;

namespace pvt = povu::types::genomics;
namespace pg = povu::genomics;
namespace pt = povu::types;
namespace bd = povu::bidirected;
namespace pgt = povu::graph_types;
namespace pst = povu::spanning_tree;
namespace pgt = povu::graph_types;
namespace pvtr = povu::tree;


namespace pic = povu::io::common;
namespace piv = povu::io::to_vcf;

namespace pcs = povu::subcommands::common;

void do_call(core::config &app_config);
}

#endif
