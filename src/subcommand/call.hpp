#ifndef PV_SUBCOMMANDS_CALL_HPP
#define PV_SUBCOMMANDS_CALL_HPP

#include <fstream>
#include <set>
#include <sstream>
#include <thread>
#include <vector>

#include "../../include/graph/bidirected.hpp"
#include "../io/to_vcf.hpp"
#include "../../include/common/types/compat.hpp"
#include "../../include/genomics/variants.hpp"
#include "../cli/app.hpp"
#include "../cli/cli.hpp"
#include "../io/pvst.hpp"
#include "./common.hpp"

namespace povu::subcommands::call {
namespace fs = std::filesystem;

namespace pvt = povu::types::genomics;
namespace pg = povu::variants;
namespace pt = povu::types;
namespace bd = povu::bidirected;
namespace pgt = povu::types::graph;
namespace pst = povu::spanning_tree;
namespace pc = povu::constants;
namespace pvtr = povu::tree;
namespace pic = povu::io::common;
namespace piv = povu::io::to_vcf;
namespace pcs = povu::subcommands::common;

void read_pvsts(const core::config &app_config, std::vector<pvtr::Tree> &pvsts);
pt::status_t get_refs(core::config &app_config);
std::vector<std::string> filter_paths_by_prefix(const core::config &app_config);
void do_call(core::config &app_config);
} // povu::subcommands::call
#endif
