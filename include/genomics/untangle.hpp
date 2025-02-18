#ifndef POVU_UNTANGLE_HPP
#define POVU_UNTANGLE_HPP

#include <cstddef>
#include <iostream>
#include <string>
#include <vector>

#include "../../src/cli/app.hpp"
#include "../common/genomics.hpp"
#include "../common/types.hpp"
#include "../graph/bidirected.hpp"
#include "../align/align.hpp"

namespace povu::untangle {
namespace pgt = povu::graph_types;
namespace pvt = povu::types::genomics;
namespace bd = povu::bidirected;
namespace pt = povu::types;
namespace pc = povu::constants;
namespace pa = povu::align;

void untangle_flb_rovs(const bd::VG &g, std::vector<pvt::RoV> rovs);
} // namespace povu::untangle

#endif // POVU_UNTANGLE_HPP
