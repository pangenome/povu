#include "./vcf.hpp"
#include <utility>
#include <vector>


namespace povu::genomics::vcf {

// move to AT class
inline bool at_match(const pvt::AW &at1, const pvt::AW &at2) {
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

pt::idx_t comp_pos(pvt::var_type_e vt, pt::idx_t step_idx) {
  std::string fn_name = std::format("[{}::{}]", MODULE, __func__);

  switch (vt) {
  case pvt::var_type_e::del:
  case pvt::var_type_e::ins:
    return step_idx - 1;
  default:
    return step_idx;
  }
}

std::vector<pvt::AW> get_alts(pt::id_t ref_id, const pvt::Exp &rws) {
  std::string fn_name = std::format("[{}::{}]", MODULE, __func__);

  std::vector<pvt::AW> alt_ats;
  const pvt::AW &ref_at = rws.get_itn(ref_id).get_at(0);

  for (auto &[alt_ref_id, alt_itn] : rws.get_ref_itns()) {

    if (alt_ref_id == ref_id) {
      continue;
    }

    const pvt::AW &alt_at = alt_itn.get_ats().front();
    //std::cerr << fn_name << " alt " << alt_at.as_str() << "\n";
    if (at_match(ref_at, alt_at)) {
      continue;
    }
    alt_ats.push_back(alt_itn.get_ats().front());
  }

  return alt_ats;
};

std::vector<pvt::AW> get_alts_tangled(const bd::VG &g, pt::id_t ref_id,
                                      pt::idx_t at_idx, const pvt::Exp &rws) {
  std::vector<pvt::AW> alt_ats;

  for (auto [up, aln] : rws.get_alns()) {
    auto [up_ref1, up_ref2] = up;

    pt::idx_t at1_count = rws.get_itn(up_ref1).at_count();
    pt::idx_t at2_count = rws.get_itn(up_ref2).at_count();

    if ((up_ref1 != ref_id && up_ref2 != ref_id) || at_idx >= std::min(at1_count, at2_count) ) {
      continue;
    }

    if (aln[at_idx] == 'D' || aln[at_idx] == 'M') {
      continue;
    }
    //std::cerr << "ref " << g.get_ref_name(up_ref1) << "\n";

    auto alt_ref_id = up_ref1 == ref_id ? up_ref2 : up_ref1;
    pvt::AW alt_at = rws.get_itn(alt_ref_id).get_at(at_idx);

    // only leave the first and last elements in alt_at

    if (aln[at_idx] == 'I') {
      alt_at.get_steps_mut().erase(alt_at.get_steps_mut().begin() + 1, alt_at.get_steps_mut().end() - 1);
    }


    alt_ats.push_back(alt_at);
  }

  return alt_ats;
};

bool is_del(const pvt::AW &aw, const pgt::walk &bounds) {
  std::string fn_name = std::format("[{}::{}]", MODULE, __func__);

  if (aw.step_count() != 2) {
    return false;
  }

  for (pt::idx_t step_idx{}; step_idx < aw.step_count(); step_idx++) {

    if (aw.get_step(step_idx).get_v_id() != bounds[step_idx].v_id ||
        aw.get_step(step_idx).get_o() != bounds[step_idx].orientation) {
      return false;
    }
  }

  return true;
}

pvt::var_type_e det_var_type(const pvst::VertexBase *pvst_vtx, const pvt::AW &ref_aw,
             const pvt::AW &alt_aw) {

  // const pvst::VertexBase *pvst_vtx = exp.get_pvst_vtx_const_ptr();
  const pvst::traversal_params_t &tp = pvst_vtx->get_traversal_params();
  pgt::id_or_t alpha = tp.start;
  pgt::id_or_t omega = tp.end;

  if ( ref_aw.step_count() < alt_aw.step_count() || is_del(ref_aw, std::vector<pgt::id_or_t>{alpha, omega})) {
    return pvt::var_type_e::del;
  }
  else if (ref_aw.step_count() > alt_aw.step_count() || is_del(alt_aw, std::vector<pgt::id_or_t>{alpha, omega})) {
    return pvt::var_type_e::ins;
  }
  else {
    return pvt::var_type_e::sub;
  }
}

void handle_allele_walk() {}

/**
 * @brief
 *
 * if not tangled itn has an AT count of 1
 */
void add_vcf_recs(const bd::VG &g, const pvt::Exp &exp, pvt::VcfRecIdx &vcf_recs) {
  std::string fn_name = std::format("[{}::{}]", MODULE, __func__);

  std::string id = exp.id();

  for (pt::id_t ref_ref_id : exp.get_ref_ids()) {
    // we know it has only one walk
    // because it is not tangled the ref_itm has only one allele walk
    pvt::Itn ref_itn = exp.get_itn(ref_ref_id);

    const pvt::AW &ref_aw = ref_itn.get_ats().front();

    const pvt::AS &s = ref_aw.get_step(1);

    std::vector<pvt::AW> alt_aws;

    std::set<pt::idx_t> alt_walks_covered;

    std::map<pvt::var_type_e, pvt::VcfRec> var_type_to_vcf_rec;

    for (pt::id_t alt_ref_id : exp.get_ref_ids()) {
      if (ref_ref_id == alt_ref_id) {
        continue;
      }

      pvt::Itn alt_itn = exp.get_itn(alt_ref_id);
      pvt::AW alt_aw = alt_itn.get_ats().front();

      // if they are the same walk idx there's no point in calling variants on them
      if (ref_aw.get_walk_idx() == alt_aw.get_walk_idx() || alt_walks_covered.contains(alt_aw.get_walk_idx())) {
        continue;
      }

      pvt::var_type_e var_typ = det_var_type(exp.get_pvst_vtx_const_ptr(), ref_aw, alt_aw);
      pt::idx_t pos = s.get_step_idx();
      //comp_pos(var_typ, s.get_step_idx());

      // inserts in-place if the key does not exist, does nothing if the key exists
      var_type_to_vcf_rec.try_emplace(var_typ, pvt::VcfRec{ref_ref_id, pos, id, ref_aw, {}, var_typ, false});

      var_type_to_vcf_rec.at(var_typ).append_alt_at(std::move(alt_aw));

      alt_walks_covered.insert(alt_aw.get_walk_idx());
    }

    for (auto &[_, r] : var_type_to_vcf_rec) {
      vcf_recs.add_rec(ref_ref_id, std::move(r));
    }
  }
}

void add_vcf_recs_tangled(const bd::VG &g, const pvt::Exp &exp, pvt::VcfRecIdx &vcf_recs) {
  std::string fn_name = std::format("[{}::{}]", MODULE, __func__);

  std::string id = exp.id();

  // std::cerr << fn_name << " RoV " << id << " ref count: " << exp.get_ref_ids().size() << "\n";

  for (pt::id_t ref_ref_id : exp.get_ref_ids()) {
    // we know it has only one walk
    // because it is not tangled the ref_itm has only one allele walk
    pvt::Itn ref_itn = exp.get_itn(ref_ref_id);

    std::vector<pvt::AW> alt_aws;

    std::set<pt::idx_t> alt_walks_covered;

    std::map<std::pair<pt::idx_t, pvt::var_type_e>, pvt::VcfRec> var_type_to_vcf_rec;

    for (pt::id_t alt_ref_id : exp.get_ref_ids()) {
      if (ref_ref_id == alt_ref_id) {
        continue;
      }

      pvt::Itn alt_itn = exp.get_itn(alt_ref_id);

      // std::cerr << fn_name << " RoV " << id << " ref " << g.get_ref_name(ref_ref_id) << " alt " << g.get_ref_name(alt_ref_id) << "\n";

      const std::string &aln = exp.get_aln(ref_ref_id, alt_ref_id);

      //std::cerr << "aln: " << aln << "\n";

      // std::cerr << "ref_itn:\n";
      // for (const auto &aw : ref_itn.get_ats()) {
      //   std::cerr << aw.as_str() << "\n";
      // }
      // std::cerr << "alt_itn:\n";
      // for (const auto &aw : alt_itn.get_ats()) {
      //   std::cerr << aw.as_str() << "\n";
      // }

      for (pt::idx_t i{}, j{}, k{}; k < aln.size() ; k++) {
        char edit_op = aln[k];
        switch (edit_op) {
        case 'M':
          i++; j++;
          continue;
        case 'D':
          j++; // skip deletion
          continue;
        case 'I':
          i++; // skip insertion
          continue;
        case 'X':
        default:
          // we only care about the positions that are 'X'
          break;
        }

        const pvt::AW &ref_aw = ref_itn.get_at(i);
        const pvt::AS &s = ref_aw.get_step(1);
        pvt::AW alt_aw = alt_itn.get_at(j);

        // if they are the same walk idx there's no point in calling variants on them
        if (ref_aw.get_walk_idx() == alt_aw.get_walk_idx() || alt_walks_covered.contains(alt_aw.get_walk_idx())) {
          continue;
        }

        // std::cerr << "ref walk idx: " << ref_aw.get_walk_idx() << " ref aw: " << ref_aw.as_str() << "\n";
        // std::cerr << "alt walk idx: " << alt_aw.get_walk_idx() << " alt aw: " << alt_aw.as_str() << "\n";

        pvt::var_type_e var_typ = det_var_type(exp.get_pvst_vtx_const_ptr(), ref_aw, alt_aw);
        pt::idx_t pos = s.get_step_idx();

        // inserts in-place if the key does not exist, does nothing if the key exists
        var_type_to_vcf_rec.try_emplace(std::make_pair(i, var_typ), pvt::VcfRec{ref_ref_id, pos, id, ref_aw, {}, var_typ, true});

        var_type_to_vcf_rec.at(std::make_pair(i, var_typ)).append_alt_at(std::move(alt_aw));

        alt_walks_covered.insert(alt_aw.get_walk_idx());

        i++; j++;
      }

    }

    std::vector<pvt::VcfRec> buf; // a buffer to hold the VCF records for this ref to sort it
    buf.reserve(var_type_to_vcf_rec.size());
    for (auto &[_, r] : var_type_to_vcf_rec) {
      buf.push_back(std::move(r));
    }
    // sort the records by position
    std::sort(buf.begin(), buf.end(), [](const pvt::VcfRec &a, const pvt::VcfRec &b) {
      return a.get_pos() < b.get_pos();
    });

    for (auto &r : buf) {
      vcf_recs.add_rec(ref_ref_id, std::move(r));
    }

    // for (auto &[_, r] : var_type_to_vcf_rec) {
    //   vcf_recs.add_rec(ref_ref_id, std::move(r));
    // }

    // if (g.get_ref_name(ref_ref_id) == "chm13__LPA__tig00000001") {
    //   exit(1);
    // }
  }

  //exit(1);
}


pvt::VcfRecIdx gen_vcf_records(const bd::VG &g, const std::vector<pvt::Exp> &exps) {
  pvt::VcfRecIdx vcf_recs;

  for (const pvt::Exp &exp : exps) {
    if (exp.is_tangled()) {
      add_vcf_recs_tangled(g, exp, vcf_recs);
    }
    else {
      add_vcf_recs(g, exp, vcf_recs);
    }

  }

  return vcf_recs;
}

}// namespace povu::genomics::vcf
