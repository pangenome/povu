#include "./vcf.hpp"

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

void add_untangled_vcf_recs(const bd::VG &g, const pvt::RefWalks &rws,
                  pvt::VcfRecIdx &vcf_recs) {


  const std::map<pt::id_t, pvt::Itn> &refs = rws.get_ref_walks();

  auto at_match = [](const pvt::AT &at1, const pvt::AT &at2) -> bool {
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

  for (const auto &[ref_id, itn] : rws.get_ref_walks()) {

    const std::vector<pvt::AT> ws = itn.get_walks();
    const pvt::AT &wf = ws.front();
    pt::idx_t pos = wf.get_steps().front().get_step_idx();

    std::string id = rws.get_flb().as_str();
    std::string format = "";

    /* determine allele traversals */

    // has only one
    const pvt::AT &ref_at = rws.get_itn(ref_id).get_walks().front();

    /* populate info field */
    std::vector<pvt::AT> alt_ats;
    for (auto &[alt_ref_id, alt_itn] : refs) {
      if (alt_ref_id == ref_id) {
        continue;
      }

      const pvt::AT &alt_at = alt_itn.get_walks().front();
      if (at_match(ref_at, alt_at)) {
        continue;
      }
      alt_ats.push_back(alt_itn.get_walks().front());
    }

    pvt::VcfRec r {ref_id, pos, id, ref_at, alt_ats, format};
    vcf_recs.add_rec(ref_id, std::move(r));
  }

}

  //std::cerr << "now has " << vcf_recs.get_recs().size() << " recs \n";



pvt::VcfRecIdx gen_vcf_records(const bd::VG &g,
                               const std::vector<pvt::RefWalks> &ref_walks) {
  pvt::VcfRecIdx vcf_recs;

  for (pt::idx_t fl_idx {}; fl_idx < ref_walks.size(); ++fl_idx) {
    const pvt::RefWalks &rw = ref_walks[fl_idx];
    if (rw.is_tangled()) {
      continue;
    }
    add_untangled_vcf_recs(g, rw, vcf_recs);
  }

  return vcf_recs;
}

}// namespace povu::genomics::vcf
