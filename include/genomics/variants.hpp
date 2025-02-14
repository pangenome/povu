#ifndef POVU_VARIANTS_HPP
#define POVU_VARIANTS_HPP

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/cli/app.hpp"
#include "../common/types.hpp"
#include "../graph/bidirected.hpp"


namespace povu::variants {
namespace pgt = povu::graph_types;
namespace bd = povu::bidirected;
namespace pt= povu::types;

void call_variants(const std::vector<pgt::flubble> &canonical_flubbles,
                   const bd::VG &bd_vg,
                   const core::config &app_config);
} // namespace povu::variants

#endif // POVU_VARIANTS_HPP
