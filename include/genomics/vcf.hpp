#ifndef POVU_GENOMICS_VCF_HPP
#define POVU_GENOMICS_VCF_HPP

#include <algorithm>
#include <vector>

#include "../common/genomics.hpp"
#include "../common/types/types.hpp"
#include "../graph/bidirected.hpp"

namespace povu::genomics::vcf {
namespace pt = povu::types;
namespace bd = povu::bidirected;
namespace pvt = povu::types::genomics;

pvt::VcfRecIdx gen_vcf_records(const bd::VG &g,
                               const std::vector<pvt::RefWalks> &ref_walks);

} // namespace povu::genomics::vcf

#endif // POVU_GENOMICS_VCF_HPP
