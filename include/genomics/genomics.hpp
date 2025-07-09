#ifndef POVU_GENOMICS_HPP
#define POVU_GENOMICS_HPP

#include "../common/types/types.hpp"
#include "../graph/bidirected.hpp"
#include "./untangle.hpp"
#include "./variants.hpp"
#include "./vcf.hpp"

namespace povu::genomics {
namespace pv = povu::variants;
namespace put = povu::untangle;
namespace pgt = povu::types::graph;
namespace bd = povu::bidirected;
namespace pvt = povu::types::genomics;
namespace pgv = povu::genomics::vcf;
namespace pt = povu::types;

pvt::VcfRecIdx gen_vcf_rec_map(const std::vector<pgt::flubble> &canonical_flubbles,
                               const bd::VG &g, const core::config &app_config);

} // namespace povu::genomics

#endif
