#include "./vcf.hpp"
#include <vector>

namespace povu::genomics::vcf {

// move to AT class
inline bool at_match(const pvt::AT &at1, const pvt::AT &at2) {
  if (at1.get_steps().size() != at2.get_steps().size()) {
    return false;
  }

  pt::idx_t step_count = at1.step_count();

  for (pt::idx_t step_idx {}; step_idx < step_count; ++step_idx) {
    if (at1.get_step(step_idx).get_v_id() != at2.get_step(step_idx).get_v_id()
        || at1.get_step(step_idx).get_o() != at2.get_step(step_idx).get_o()) {
      return false;
    }
  }

  return true;
};

std::vector<pvt::AT>  get_alts (pt::id_t ref_id, const pvt::RefWalks &rws) {
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

std::vector<pvt::AT> get_alts_tangled(const bd::VG &g, pt::id_t ref_id,
                                      pt::idx_t at_idx,
                                      const pvt::RefWalks &rws) {
  std::vector<pvt::AT> alt_ats;

  for (auto [up, aln] : rws.get_alns()) {
    auto [up_ref1, up_ref2] = up;

    auto alt_id = up_ref1 == ref_id ? up_ref2 : up_ref1;

    if (aln.size() <= at_idx || aln[at_idx] == 'M') {
      continue;
    }

    // TODO: call from rws
    pvt::It alt_itn = rws.get_itn(alt_id);
    pvt::AT alt_at = alt_itn.get_at(at_idx);

    alt_ats.push_back(alt_at);
  }

  return alt_ats;
};

/**
 * @brief
 *
 * if not tangled itn has an AT count of 1
 */
void add_vcf_recs(const bd::VG &g, const pvt::RefWalks &rws,
                  pvt::VcfRecIdx &vcf_recs) {

  std::string id = rws.get_flb().as_str();

  bool dbg = false;
  if (dbg && rws.is_tangled()) {
    std::cerr << "flb: " << id<< "\n";
  }

  for (const auto &[ref_id, itn] : rws.get_ref_itns()) {

    for (pt::idx_t at_idx {}; at_idx < itn.at_count(); at_idx++) {

      const pvt::AT &ref_at = itn.get_at(at_idx);
      const pvt::Step &s = ref_at.get_step(1);
      pt::id_t pos = ref_at.is_del() ? s.get_step_idx() - 1 : s.get_step_idx();

      std::vector<pvt::AT> alt_ats = rws.is_tangled() ?
        get_alts_tangled(g, ref_id, at_idx, rws) : get_alts(ref_id, rws);

      if (dbg && rws.is_tangled()) {
        std::cerr << std::format("ref: {} \nAlts: ", ref_at.as_str());
        for (pt::idx_t at_idx{}; at_idx < alt_ats.size(); at_idx++) {
          const pvt::AT &at = alt_ats[at_idx];
          std::cerr << at.as_str();
          if (at_idx < itn.at_count() - 1) {
            std::cerr << ", ";
          }
        }
        std::cerr << "\n";
      }

      pvt::VcfRec r{ref_id, pos, id, ref_at, alt_ats};
      vcf_recs.add_rec(ref_id, std::move(r));
    }
  }

  if (dbg && rws.is_tangled()) {
    for (auto &[k, v] : rws.get_alns()) {
      std::cerr << std::format("[{},{}] aln {}\n", g.get_ref_name(k.l), g.get_ref_name(k.r), v);
    }
    exit(1);
  }
}

pvt::VcfRecIdx gen_vcf_records(const bd::VG &g, const std::vector<pvt::RefWalks> &ref_walks) {
  pvt::VcfRecIdx vcf_recs;

  for (pt::idx_t fl_idx {}; fl_idx < ref_walks.size(); ++fl_idx) {
    // a refWalk applies to a single flubble
    const pvt::RefWalks &rw = ref_walks[fl_idx];
    add_vcf_recs(g, rw, vcf_recs);
  }

  return vcf_recs;
}

}// namespace povu::genomics::vcf
