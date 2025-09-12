#include "./vcf.hpp"


namespace povu::genomics::vcf {



var_type_e det_var_type(const pga::allele_slice_t &ref_allele_slice,
                        const pga::allele_slice_t &alt_allele_slice) {
  if (ref_allele_slice.len < alt_allele_slice.len) {
    return var_type_e::del;
  }
  else if (ref_allele_slice.len > alt_allele_slice.len) {
    return var_type_e::ins;
  }
  else {
    return var_type_e::sub;
  }
}

pt::idx_t comp_pos(const pga::allele_slice_t &ref_allele_slice, var_type_e variant_type) {
  const pgt::ref_walk_t rw = *ref_allele_slice.ref_walk;
  pt::idx_t locus = rw[ref_allele_slice.ref_start_idx + 1].locus;

  switch (variant_type) {
  case var_type_e::del:
  case var_type_e::ins:
    return locus - 1;
  default:
    return locus;
  }
}

std::vector<pt::op_t<pt::idx_t>> get_call_itn_idxs(const pga::Exp &exp, pt::id_t ref_w_ref_id, pt::id_t alt_w_ref_id) {
  std::vector<pt::op_t<pt::idx_t>> idxs_to_call;

  if (!exp.is_tangled() || !exp.has_aln(ref_w_ref_id, alt_w_ref_id)) {
    idxs_to_call.emplace_back(0, 0);
    return idxs_to_call;
  }

  const std::string &aln = exp.get_aln(ref_w_ref_id, alt_w_ref_id);

  for (pt::idx_t i{}, j{}, k{}; k < aln.size(); k++) {
    char edit_op = aln[k];
    switch (edit_op) {
    case 'M':
      i++;
      j++;
      break;
    case 'D':
      j++; // skip deletion
      break;
    case 'I':
      i++; // skip insertion
      break;
    case 'X':
    default:
      // we only care about the positions that are 'X'
      idxs_to_call.emplace_back(i, j);
      i++;
      j++;
      break;
    }
  }

  return idxs_to_call;
}

// ref id to a list of vcf records for that ref
std::map<pt::idx_t, std::vector<VcfRec>>
gen_exp_vcf_recs(const bd::VG &g, const pga::Exp &exp, const std::set<pt::id_t> &to_call_ref_ids) {
  std::map<pt::idx_t, std::vector<VcfRec>> exp_vcf_recs;

  if (exp.get_rov() == nullptr) {
    ERR("RoV pointer is null");
    std::exit(EXIT_FAILURE);
  }
  const pgg::RoV &rov = *(exp.get_rov());

  if (exp.get_pvst_vtx_const_ptr() == nullptr) {
    ERR("pvst vertex pointer is null");
    std::exit(EXIT_FAILURE);
  }
  const pvst::VertexBase *pvst_vtx_ptr = exp.get_pvst_vtx_const_ptr();

  auto pre_comp_ref_pairs = [&]() -> std::vector<pt::op_t<pt::id_t>> {
    std::vector<pt::op_t<pt::id_t>> ref_pairs;
    std::set<pt::id_t> exp_ref_ids = exp.get_ref_ids();
    for (pt::id_t ref_ref_id : exp_ref_ids) {
      if (!pv_cmp::contains(to_call_ref_ids, ref_ref_id)) {
        continue;
      }
      for (pt::id_t alt_ref_id : exp_ref_ids) {
        if (ref_ref_id != alt_ref_id) {
          ref_pairs.emplace_back(ref_ref_id, alt_ref_id);
        }
      }
    }
    return ref_pairs;
  };

  std::map<std::pair<pt::idx_t, var_type_e>, VcfRec> var_type_to_vcf_rec;

  for (auto [ref_ref_id, alt_ref_id] : pre_comp_ref_pairs()) {
    const pga::itn_t &ref_itn = exp.get_itn(ref_ref_id);
    const pga::itn_t &alt_itn = exp.get_itn(alt_ref_id);

    for (auto [i, j] : get_call_itn_idxs(exp, ref_ref_id, alt_ref_id)) {
      const pga::allele_slice_t &ref_allele_slice = ref_itn.get_at(i);
      const pga::allele_slice_t &alt_allele_slice = alt_itn.get_at(j);

      pt::idx_t ref_walk_idx = ref_allele_slice.walk_idx;
      pt::idx_t alt_walk_idx = alt_allele_slice.walk_idx;

      // TODO: check if start and len of the walks as well
      if (ref_walk_idx == alt_walk_idx) { // this means they are from the same walk, skip
        continue;
      }

      pt::idx_t ref_walk_ref_count = exp.get_ref_idxs_for_walk(ref_walk_idx).size();
      pt::idx_t alt_walk_ref_count = exp.get_ref_idxs_for_walk(alt_walk_idx).size();

      var_type_e variant_type = det_var_type(ref_allele_slice, alt_allele_slice);

      std::pair<pt::idx_t, var_type_e> key = std::make_pair(ref_ref_id, variant_type);

      // if it does not exist create a variant type for it and add to var_type_to_vcf_rec
      if (!pv_cmp::contains(var_type_to_vcf_rec, key)) {
        pt::idx_t pos = comp_pos(ref_allele_slice, variant_type);
        VcfRec vcf_rec {
          ref_ref_id,
          pos,
          exp.id(),
          ref_allele_slice,
          pvst_vtx_ptr->get_height(),
          variant_type,
          exp.is_tangled(),
          ref_walk_ref_count,
          g.get_blank_genotype_cols()};

        var_type_to_vcf_rec.emplace(key, std::move(vcf_rec));
      }

      VcfRec &curr_vcf_rec = var_type_to_vcf_rec.at(key);
      const pt::idx_t alt_allele_col_idx = curr_vcf_rec.append_alt_at(alt_allele_slice, alt_walk_ref_count);
    }
  }

  for (auto &[k, r] : var_type_to_vcf_rec) {
    auto [ref_ref_id, _] = k;
    exp_vcf_recs[ref_ref_id].emplace_back(std::move(r));
  }

  return exp_vcf_recs;
}

VcfRecIdx gen_vcf_records(const bd::VG &g, const std::vector<pga::Exp> &exps,
                          const std::set<pt::id_t> &to_call_ref_ids) {
  VcfRecIdx vcf_recs;

  for (pt::idx_t i {}; i < exps.size(); ++i) {
    const pga::Exp &exp = exps[i];
    if (exp.get_pvst_vtx_const_ptr() == nullptr) {
      ERR("pvst vertex pointer is null {}", exp.get_rov()->as_str());
      std::exit(EXIT_FAILURE);
    }
    std::map<pt::idx_t, std::vector<VcfRec>> exp_vcf_recs = gen_exp_vcf_recs(g, exp, to_call_ref_ids);
    vcf_recs.ensure_append_recs(std::move(exp_vcf_recs));
  }

  return vcf_recs;
}

}// namespace povu::genomics::vcf
