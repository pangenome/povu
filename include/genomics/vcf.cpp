#include "./vcf.hpp"
#include <string>
#include <vector>

namespace povu::genomics::vcf {

void gen_flb_vcf_recs(pvt::VcfRecIdx &vcf_recs , const pvt::RefWalks &ref_walks) {

  for (const auto &[ref_id, itn] : ref_walks.get_ref_walks()) {
    //std::string chrom = g.get_ref_name(ref_id);

    const std::vector<pvt::Walk> ws = itn.get_walks();
    const pvt::Walk &wf = ws.front();
    pt::idx_t pos = wf.get_steps().front().get_step_idx();
    std::string id = "";
    std::string format = "";

    for (const pvt::Walk &ref_at: ws) {
      std::vector<pvt::Walk> alt_ats;
      vcf_recs.add_rec(ref_id,
                       pvt::VcfRec{ref_id, pos, id, ref_at, alt_ats, format});
    }
  }
}

pvt::VcfRecIdx to_vcf_records(const std::vector<pvt::RefWalks> &ref_walks) {

  pvt::VcfRecIdx vcf_recs;
  for (const auto &ref_walk : ref_walks) {
    gen_flb_vcf_recs(vcf_recs, ref_walk);
  }

  return vcf_recs;
}

}// namespace povu::genomics::vcf
