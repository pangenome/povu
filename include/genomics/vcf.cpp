#include "./vcf.hpp"
#include <string>

namespace povu::genomics::vcf {

void gen_vcfs(const bd::VG &g, const std::vector<pvt::RefWalks> &ref_walks) {

  for (const auto &[ref_id, ref_walk] : ref_walks.get_ref_walks()) {
    std::string chrom = g.get_ref_name(ref_id);

  }
}

}// namespace povu::genomics::vcf
