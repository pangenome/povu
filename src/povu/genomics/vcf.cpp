#include "povu/genomics/vcf.hpp"

#include <iostream>
#include <iterator> // for pair
#include <vector>

#include "povu/common/core.hpp"
#include "povu/common/utils.hpp"
#include "povu/genomics/allele.hpp" // for Exp, allele_slice_t, itn_t
#include "povu/genomics/graph.hpp"  // for RoV
#include "povu/graph/pvst.hpp"	    // for VertexBase

namespace povu::genomics::vcf
{
namespace pvst = povu::pvst;
namespace pgg = povu::genomics::graph;

constexpr pvr::var_type_e ins = pvr::var_type_e::ins;
constexpr pvr::var_type_e del = pvr::var_type_e::del;
constexpr pvr::var_type_e sub = pvr::var_type_e::sub;

pvr::var_type_e det_var_type(const pga::allele_slice_t &ref_allele_slice,
			     const pga::allele_slice_t &alt_allele_slice)
{
	if (ref_allele_slice.len < alt_allele_slice.len) {
		return ins;
	}
	else if (ref_allele_slice.len > alt_allele_slice.len) {
		return del;
	}
	else {
		return sub;
	}
}

pt::idx_t comp_pos(const pga::allele_slice_t &ref_allele_slice,
		   pvr::var_type_e variant_type)
{
	// const pgt::ref_walk_t rw = *ref_allele_slice.ref_walk;
	// pt::idx_t locus = rw[ref_allele_slice.ref_start_idx + 1].locus;

	pt::idx_t locus =
		ref_allele_slice.get_locus(ref_allele_slice.ref_start_idx + 1);

	switch (variant_type) {
	case del:
	case ins:
		return locus - 1;
	default:
		return locus;
	}
}

std::vector<pt::op_t<pt::idx_t>> get_call_itn_idxs(const pga::Exp &exp,
						   pt::id_t ref_w_ref_id,
						   pt::id_t alt_w_ref_id)
{
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

std::vector<pt::op_t<pt::u32>>
pre_comp_ref_pairs(const pga::Exp &exp, const std::set<pt::u32> &call_ref_ids)
{
	std::vector<pt::op_t<pt::id_t>> ref_pairs;
	std::set<pt::u32> exp_ref_ids = exp.get_ref_ids();

	for (pt::id_t ref_ref_id : exp_ref_ids) {
		if (!pv_cmp::contains(call_ref_ids, ref_ref_id))
			continue;

		for (pt::id_t alt_ref_id : exp_ref_ids)
			if (ref_ref_id != alt_ref_id)
				ref_pairs.emplace_back(ref_ref_id, alt_ref_id);
	}

	return ref_pairs;
}

bool non_varying(const pga::allele_slice_t &ref_allele_slice,
		 const pga::allele_slice_t &alt_allele_slice)
{
	pt::idx_t ref_walk_idx = ref_allele_slice.walk_idx;
	pt::idx_t alt_walk_idx = alt_allele_slice.walk_idx;

	auto ref_sl_or = ref_allele_slice.slice_or;
	auto alt_sl_or = alt_allele_slice.slice_or;

	return (ref_walk_idx == alt_walk_idx && ref_sl_or == alt_sl_or);
}

// ref id to a list of vcf records for that ref
std::map<pt::idx_t, std::vector<VcfRec>>
gen_exp_vcf_recs(const bd::VG &g, const pga::Exp &exp,
		 const std::set<pt::id_t> &to_call_ref_ids)
{
#ifdef DEBUG
	if (exp.get_rov() == nullptr) {
		ERR("RoV pointer is null");
		std::exit(EXIT_FAILURE);
	}

	if (exp.get_pvst_vtx_const_ptr() == nullptr) {
		ERR("pvst vertex pointer is null");
		std::exit(EXIT_FAILURE);
	}
#endif

	std::map<pt::idx_t, std::vector<VcfRec>> exp_vcf_recs;
	const pvr::RoV &rov = *(exp.get_rov());
	const pvst::VertexBase *pvst_vtx_ptr = exp.get_pvst_vtx_const_ptr();

	std::map<std::tuple<pt::idx_t, pt::u32, pvr::var_type_e>, VcfRec>
		var_type_to_vcf_rec;

	std::tuple<pt::idx_t, pt::u32, pvr::var_type_e> key;

	std::vector<pt::op_t<pt::id_t>> ref_pairs =
		pre_comp_ref_pairs(exp, to_call_ref_ids);
	for (auto [ref_ref_id, alt_ref_id] : ref_pairs) {

		const pga::itn_t &ref_itn = exp.get_itn(ref_ref_id);
		const pga::itn_t &alt_itn = exp.get_itn(alt_ref_id);

		if (!exp.has_ref(ref_ref_id))
			continue;

		for (auto [i, j] :
		     get_call_itn_idxs(exp, ref_ref_id, alt_ref_id)) {
			pga::allele_slice_t ref_allele_slice =
				ref_itn.get_at(i);
			pga::allele_slice_t alt_allele_slice =
				alt_itn.get_at(j);

			pt::idx_t ref_walk_idx = ref_allele_slice.walk_idx;
			pt::idx_t alt_walk_idx = alt_allele_slice.walk_idx;

			if (non_varying(ref_allele_slice, alt_allele_slice))
				continue;

			if (ref_allele_slice == alt_allele_slice)
				continue;
			// TODO: check if start and len of the walks as
			// well this means they are from the same walk,
			// skip
			if (ref_walk_idx == alt_walk_idx)
				continue;

			pt::idx_t ref_walk_ref_count =
				exp.get_ref_idxs_for_walk(ref_walk_idx).size();
			pt::idx_t alt_walk_ref_count =
				exp.get_ref_idxs_for_walk(alt_walk_idx).size();

			pvr::var_type_e variant_type = det_var_type(
				ref_allele_slice, alt_allele_slice);

			key = {ref_ref_id, ref_walk_idx, variant_type};

			// if it does not exist create a variant type for it and
			// add to var_type_to_vcf_rec
			if (!pv_cmp::contains(var_type_to_vcf_rec, key)) {

				pt::idx_t pos = comp_pos(ref_allele_slice,
							 variant_type);

				VcfRec vcf_rec{ref_ref_id,
					       pos,
					       exp.id(),
					       ref_allele_slice,
					       pvst_vtx_ptr->get_height(),
					       variant_type,
					       exp.is_tangled(),
					       ref_walk_ref_count,
					       g.get_blank_genotype_cols()};

				var_type_to_vcf_rec.emplace(key,
							    std::move(vcf_rec));
			}

			VcfRec &curr_vcf_rec = var_type_to_vcf_rec.at(key);
			const pt::idx_t alt_allele_col_idx =
				curr_vcf_rec.append_alt_at(alt_allele_slice,
							   alt_walk_ref_count);
		}
	}

	for (auto &[k, r] : var_type_to_vcf_rec) {
		auto [ref_ref_id, _, __] = k;
		exp_vcf_recs[ref_ref_id].emplace_back(std::move(r));
	}

	return exp_vcf_recs;
}

VcfRecIdx gen_vcf_records(const bd::VG &g, const std::vector<pga::Exp> &exps,
			  const std::set<pt::id_t> &to_call_ref_ids)
{
	VcfRecIdx vcf_recs;

	for (pt::idx_t i{}; i < exps.size(); ++i) {
		const pga::Exp &exp = exps[i];
		if (exp.get_pvst_vtx_const_ptr() == nullptr) {
			ERR("pvst vertex pointer is null {}",
			    exp.get_rov()->as_str());
			std::exit(EXIT_FAILURE);
		}

		std::map<pt::idx_t, std::vector<VcfRec>> exp_vcf_recs =
			gen_exp_vcf_recs(g, exp, to_call_ref_ids);
		vcf_recs.ensure_append_recs(std::move(exp_vcf_recs));
	}

	return vcf_recs;
}

} // namespace povu::genomics::vcf
