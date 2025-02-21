#include "./vcf.hpp"
#include <string>
#include <utility>
#include <vector>

namespace povu::genomics::vcf {



std::vector<pvt::AT> get_alts(pt::id_t ref_id, pt::idx_t at_idx, const pvt::RefWalks &rws) {
  std::vector<pvt::AT> alt_ats;
  for (auto [up, aln] : rws.get_alns()) {
    auto [up_ref1, up_ref2] = up;
    if (up_ref1 != ref_id || up_ref2 != ref_id || aln.size() <= at_idx || aln[at_idx] == 'M') {
      continue;
    }

    auto alt_id = up_ref1 == ref_id ? up_ref2 : up_ref1;
    pvt::It alt_itn = rws.get_itn(alt_id);
    pvt::AT alt_at = alt_itn.get_walk_by_step_idx(at_idx);
    alt_ats.push_back(alt_at);
  }

  return alt_ats;
}

void add_vcf_recs(const pvt::RefWalks &rws, pvt::VcfRecIdx &vcf_recs) {
  for (const auto &[ref_id, itn] : rws.get_ref_walks()) {
    const std::vector<pvt::AT> ws = itn.get_walks();
    const pvt::AT &wf = ws.front();
    pt::idx_t pos = wf.get_steps().front().get_step_idx();
    std::string id = "";
    std::string format = "";

    for (pt::idx_t at_idx {}; at_idx < ws.size(); ++at_idx) {
      std::vector<pvt::AT> alt_ats = get_alts(ref_id, at_idx, rws);
      if (alt_ats.empty()) {
        continue;
      }

      const pvt::AT &ref_at = ws[at_idx];
      pvt::VcfRec r {ref_id, pos, id, ref_at, alt_ats, format};
      vcf_recs.add_rec(ref_id, std::move(r));
    }
  }
}

pvt::VcfRecIdx gen_vcf_records(const std::vector<pvt::RefWalks> &ref_walks) {

  pvt::VcfRecIdx vcf_recs;

  for (pt::idx_t fl_idx {}; fl_idx < ref_walks.size(); ++fl_idx) {
    const pvt::RefWalks &rw = ref_walks[fl_idx];
    add_vcf_recs(rw, vcf_recs);
  }

  return vcf_recs;
}

}// namespace povu::genomics::vcf
