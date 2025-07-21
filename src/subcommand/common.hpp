#ifndef PV_COMMON_SUBCOMMANDS_HPP
#define PV_COMMON_SUBCOMMANDS_HPP

#include "../../include/common/types/types.hpp"
#include "../../include/graph/bidirected.hpp"
#include "../cli/app.hpp"
#include "../io/from_gfa.hpp"

namespace povu::subcommands::common {
namespace bd = povu::bidirected;
namespace pt = povu::types;

/* ------ common (or utility) functions ------- */
bd::VG *get_vg(const core::config &app_config);
}

#endif
