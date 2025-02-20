#ifndef POVU_GENOMICS_HPP
#define POVU_GENOMICS_HPP

#include "../common/types.hpp"
#include "../graph/bidirected.hpp"
#include "./untangle.hpp"
#include "./variants.hpp"
#include "./vcf.hpp"

namespace povu::genomics {
namespace pv = povu::variants;
namespace put = povu::untangle;
namespace pgt = povu::graph_types;
namespace bd = povu::bidirected;
namespace pvt = povu::types::genomics;
namespace pgv = povu::genomics::vcf;

void to_vcf(const std::vector<pgt::flubble> &canonical_flubbles,
            const bd::VG &g,
            const core::config &app_config);

} // namespace povu::genomics

#endif
