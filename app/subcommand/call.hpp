#ifndef PV_SUBCOMMANDS_CALL_HPP
#define PV_SUBCOMMANDS_CALL_HPP

#include <fstream>
#include <set>
#include <sstream>
#include <string_view>
#include <thread>
#include <vector>


#include "../../include/graph/spanning_tree.hpp"
#include "../../include/graph/bidirected.hpp"
#include "../../include/io/to_vcf.hpp"
#include "../../include/common/compat.hpp"
#include "../../include/genomics/genomics.hpp"
#include "../cli/app.hpp"
#include "../cli/cli.hpp"
#include "../../include/io/from_pvst.hpp"
#include "./common.hpp"

namespace povu::subcommands::call {

constexpr std::string_view MODULE = "povu::subcommands::call";

namespace fs = std::filesystem;

namespace pvt = povu::types::genomics;
namespace pg = povu::genomics;
namespace bd = povu::bidirected;
namespace pgt = povu::types::graph;
namespace pst = povu::spanning_tree;
namespace pc = povu::constants;
namespace pvtr = povu::tree;
namespace pic = povu::io::common;
namespace piv = povu::io::to_vcf;


void read_pvsts(const core::config &app_config, std::vector<pvtr::Tree> &pvsts);
pt::status_t get_refs(core::config &app_config);
std::vector<std::string> filter_paths_by_prefix(const core::config &app_config);
void do_call(core::config &app_config);
} // povu::subcommands::call
#endif
