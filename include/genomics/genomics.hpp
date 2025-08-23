#ifndef POVU_GENOMICS_HPP
#define POVU_GENOMICS_HPP

#include "../common/types/genomics.hpp"
#include "../graph/bidirected.hpp"
#include "../../include/graph/tree.hpp"
#include "../common/compat.hpp"
#include "./allele.hpp"
#include "./graph.hpp"
#include "./untangle.hpp"
#include "./vcf.hpp"
#include "../common/utils.hpp"
#include "../common/log.hpp"

namespace povu::genomics {
inline constexpr std::string_view MODULE = "povu::variants";

namespace put = povu::genomics::untangle;
namespace pgt = povu::types::graph;
namespace pvt = povu::types::genomics;
namespace pgv = povu::genomics::vcf;
namespace pvtr = povu::tree;
namespace pvst = povu::types::pvst;

pvt::VcfRecIdx gen_vcf_rec_map(const std::vector<pvtr::Tree> &pvsts, bd::VG &g);

} // namespace povu::genomics

#endif
