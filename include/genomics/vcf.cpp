#include "./vcf.hpp"
#include <vector>

namespace povu::genomics::vcf {

// move to AT class
inline bool at_match(const pvt::AT &at1, const pvt::AT &at2) {
  if (at1.get_steps().size() != at2.get_steps().size()) {
    return false;
  }

  for (pt::idx_t i {}; i < at1.get_steps().size(); ++i) {
    if (at1.get_steps()[i].get_v_id() != at2.get_steps()[i].get_v_id() ||
        at1.get_steps()[i].get_o() != at2.get_steps()[i].get_o()) {
      return false;
    }
  }

  return true;
};

/**
 * @brief
 *
 * if not tangled itn has an AT count of 1
 */
void add_vcf_recs(const pvt::RefWalks &rws, pvt::VcfRecIdx &vcf_recs) {

  auto get_alts = [](pt::id_t ref_id, const pvt::RefWalks &rws) {
    std::vector<pvt::AT> alt_ats;

    const pvt::AT &ref_at = rws.get_itn(ref_id).get_at(0);
    for (auto &[alt_ref_id, alt_itn] : rws.get_ref_itns()) {
      if (alt_ref_id == ref_id) {
        continue;
      }

      const pvt::AT &alt_at = alt_itn.get_ats().front();
      if (at_match(ref_at, alt_at)) {
        continue;
      }
      alt_ats.push_back(alt_itn.get_ats().front());
    }

    return alt_ats;
  };

  auto get_alts_tangled = [](pt::id_t ref_id, pt::idx_t at_idx, const pvt::RefWalks &rws) -> std::vector<pvt::AT> {
    std::vector<pvt::AT> alt_ats;

    for (auto [up, aln] : rws.get_alns()) {
      auto [up_ref1, up_ref2] = up;
      if (up_ref1 != ref_id || up_ref2 != ref_id || aln.size() <= at_idx || aln[at_idx] == 'M') {
        continue;
      }

      auto alt_id = up_ref1 == ref_id ? up_ref2 : up_ref1;
      // TODO: call from rws
      pvt::It alt_itn = rws.get_itn(alt_id);
      pvt::AT alt_at = alt_itn.get_at(at_idx);
      alt_ats.push_back(alt_at);
    }

    return alt_ats;
  };

  for (const auto &[ref_id, itn] : rws.get_ref_itns()) {

    for (pt::idx_t at_idx{}; at_idx < itn.at_count(); ++at_idx) {

      const pvt::AT &at = itn.get_at(at_idx);
      const pvt::Step &s = at.get_step(1);
      pt::id_t pos = at.is_del() ? s.get_step_idx() - 1 : s.get_step_idx();

      std::string id = rws.get_flb().as_str();
      //std::string format = ""; // TODO: not needed here

      const pvt::AT &ref_at = at;
      std::vector<pvt::AT> alt_ats =
        rws.is_tangled() ?  get_alts(ref_id, rws) : get_alts_tangled(ref_id, at_idx, rws);

      pvt::VcfRec r{ref_id, pos, id, ref_at, alt_ats};
      vcf_recs.add_rec(ref_id, std::move(r));
    }
  }
}

pvt::VcfRecIdx gen_vcf_records(const bd::VG &g, const std::vector<pvt::RefWalks> &ref_walks) {
  pvt::VcfRecIdx vcf_recs;

  for (pt::idx_t fl_idx {}; fl_idx < ref_walks.size(); ++fl_idx) {
    // a refWalk applies to a single flubble
    const pvt::RefWalks &rw = ref_walks[fl_idx];
    add_vcf_recs(rw, vcf_recs);
  }

  return vcf_recs;
}

}// namespace povu::genomics::vcf
