#ifndef POVU_VARIANTS_HPP
#define POVU_VARIANTS_HPP

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "../../include/graph/tree.hpp"
#include "../../src/cli/app.hpp"
#include "../common/types/genomics.hpp"
#include "../common/types/types.hpp"
#include "../graph/bidirected.hpp"
#include "../common/types/pvst.hpp"
#include "../common/graph_utils.hpp"

namespace povu::variants {
inline constexpr std::string_view MODULE = "povu::variants";

namespace pgt = povu::types::graph;
namespace bd = povu::bidirected;
namespace pt= povu::types;
namespace pvt = povu::types::genomics;
namespace pvtr = povu::tree;
namespace pvst = povu::types::pvst;


std::vector<pvt::RoV> gen_rov(const std::vector<pvtr::Tree> &pvsts,
                              const bd::VG &g, const core::config &app_config);
} // namespace povu::variants


#endif // POVU_VARIANTS_HPP
