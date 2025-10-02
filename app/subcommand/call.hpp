#ifndef PV_SUBCOMMANDS_CALL_HPP
#define PV_SUBCOMMANDS_CALL_HPP

#include <cassert>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <ostream>
#include <set>
#include <sstream>
#include <string_view>
#include <thread>
#include <vector>

#include "../../include/common/bounded_queue.hpp"
#include "../../include/common/compat.hpp"
#include "../../include/common/log.hpp"
#include "../../include/common/progress.hpp"
#include "../../include/genomics/allele.hpp"
#include "../../include/genomics/genomics.hpp"
#include "../../include/genomics/vcf.hpp"
#include "../../include/graph/bidirected.hpp"
#include "../../include/graph/spanning_tree.hpp"
#include "../../include/io/from_gfa.hpp"
#include "../../include/io/from_pvst.hpp"
#include "../../include/io/to_vcf.hpp"
#include "../cli/app.hpp"
#include "../cli/cli.hpp"

namespace povu::subcommands::call
{

constexpr std::string_view MODULE = "povu::subcommands::call";

using namespace povu::progress;

namespace fs = std::filesystem;

namespace pga = povu::genomics::allele;
namespace pgv = povu::genomics::vcf;
namespace pg = povu::genomics;
namespace bd = povu::bidirected;
namespace pgt = povu::types::graph;
namespace pst = povu::spanning_tree;
namespace pc = povu::constants;
namespace pvst = povu::pvst;
namespace pic = povu::io::common;
namespace piv = povu::io::to_vcf;

void do_call(core::config &app_config);
} // namespace povu::subcommands::call
#endif
