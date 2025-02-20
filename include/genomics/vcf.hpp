#ifndef POVU_GENOMICS_VCF_HPP
#define POVU_GENOMICS_VCF_HPP

#include "../common/genomics.hpp"
#include "../graph/bidirected.hpp"

namespace povu::genomics::vcf {
namespace pvt = povu::types::genomics;
namespace bd = povu::bidirected;

void gen_vcfs(const bd::VG &g, const std::vector<pvt::RefWalks> &ref_walks);

} // namespace povu::genomics::vcf

#endif // POVU_GENOMICS_VCF_HPP
