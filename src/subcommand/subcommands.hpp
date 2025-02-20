#ifndef POVU_SUBCOMMANDS_HPP
#define POVU_SUBCOMMANDS_HPP

#include <thread> // for deconstruct

#include "../../include/algorithms/algorithms.hpp"
#include "../../include/common/types.hpp"
#include "../../include/graph/flubble_tree.hpp"
#include "../../include/graph/spanning_tree.hpp"
#include "../../include/graph/bidirected.hpp"
#include "../io/to_vcf.hpp"
#include "../cli/app.hpp"
#include "../cli/cli.hpp"
#include "../../include/genomics/genomics.hpp"
#include "../../include/common/genomics.hpp"

namespace povu::subcommands {
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

void do_call(core::config &app_config);
void do_info(const core::config &app_config);
void do_deconstruct(const core::config &app_config);

/* ------ common (or utility) functions ------- */
bd::VG *get_vg(const core::config &app_config);
}

#endif
