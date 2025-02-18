#ifndef POVU_VARIANTS_HPP
#define POVU_VARIANTS_HPP

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/cli/app.hpp"
#include "../common/types.hpp"
#include "../common/genomics.hpp"
#include "../graph/bidirected.hpp"


namespace povu::variants {
namespace pgt = povu::graph_types;
namespace bd = povu::bidirected;
namespace pt= povu::types;
namespace pvt = povu::types::genomics;

// Maximum number of steps to take from flubble start to end
const pt::idx_t MAX_FLUBBLE_STEPS {20};

std::vector<pvt::RoV> gen_rov(const std::vector<pgt::flubble> &canonical_flubbles,
                              const bd::VG &g, const core::config &app_config);
} // namespace povu::variants


#endif // POVU_VARIANTS_HPP
