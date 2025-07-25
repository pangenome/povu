#ifndef POVU_GENOMICS_HPP
#define POVU_GENOMICS_HPP

#include "../../include/graph/tree.hpp"
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
namespace pvtr = povu::tree;


pvt::VcfRecIdx gen_vcf_rec_map(const std::vector<pvtr::Tree> &pvsts,
                               const bd::VG &g, const core::config &app_config);

} // namespace povu::genomics

#endif
