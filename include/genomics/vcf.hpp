#ifndef POVU_GENOMICS_VCF_HPP
#define POVU_GENOMICS_VCF_HPP

#include <algorithm>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "../common/types/genomics.hpp"
#include "../common/compat.hpp"
#include "../common/types/graph.hpp"
#include "../graph/bidirected.hpp"

namespace povu::genomics::vcf {
inline constexpr std::string_view MODULE = "povu::genomics::vcf";

namespace bd = povu::bidirected;
namespace pvt = povu::types::genomics;
namespace pgt = povu::types::graph;
namespace pvst = povu::types::pvst;


pvt::VcfRecIdx gen_vcf_records(const bd::VG &g, const std::vector<pvt::Exp> &ref_walks);

} // namespace povu::genomics::vcf

#endif // POVU_GENOMICS_VCF_HPP
