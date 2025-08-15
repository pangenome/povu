#ifndef POVU_GENOMICS_HPP
#define POVU_GENOMICS_HPP

#include "../../include/graph/tree.hpp"
#include "./utils.hpp"
#include "../common/compat.hpp"
#include "../common/types/genomics.hpp"
#include "../graph/bidirected.hpp"
#include "./untangle.hpp"
#include "./vcf.hpp"

namespace povu::variants {
namespace put = povu::untangle;
namespace pgt = povu::types::graph;
namespace pgu = povu::genomics::utils;
namespace bd = povu::bidirected;
namespace pvt = povu::types::genomics;
namespace pgv = povu::genomics::vcf;
namespace pt = povu::types;
namespace pvtr = povu::tree;
namespace pvst = povu::types::pvst;

inline constexpr std::string_view MODULE = "povu::variants";

pvt::VcfRecIdx gen_vcf_rec_map(const std::vector<pvtr::Tree> &pvsts, const bd::VG &g);

} // namespace povu::genomics

#endif
