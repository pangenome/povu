#ifndef PV_COMMON_SUBCOMMANDS_HPP
#define PV_COMMON_SUBCOMMANDS_HPP

#include "../../include/common/compat.hpp"
#include "../../include/graph/bidirected.hpp"
#include "../cli/app.hpp"
#include "../../include/io/from_gfa.hpp"

namespace povu::subcommands::common {
namespace bd = povu::bidirected;

/* ------ common (or utility) functions ------- */
bd::VG *get_vg(const core::config &app_config);
}

#endif
