#include "povu/genomics/vcf.hpp"

#include <optional>
#include <string>
#include <vector>

#include <liteseq/refs.h> // for ref_walk, ref

#include "povu/common/core.hpp"
#include "povu/genomics/allele.hpp" // for Exp, allele_slice_t, itn_t
#include "povu/graph/pvst.hpp"	    // for VertexBase
#include "povu/tree/slice_tree.hpp" // for poi
#include "povu/variation/rov.hpp"

namespace povu::genomics::vcf
{
namespace lq = liteseq;
namespace pvst = povu::pvst;

std::pair<pt::u32, std::vector<std::vector<std::string>>>
gen_gt_data(const bd::VG &g, const std::set<pt::u32> &ref_haps,
	    const std::set<pt::u32> &alt_h_idxs)
{
	std::vector<std::vector<std::string>> gt_cols =
		g.get_blank_genotype_cols();

	for (std::vector<std::string> &col : gt_cols)
		for (std::string &gt : col)
			gt = ".";

	// number of samples with data. cols that that are non-zero
	std::set<pt::u32> ns_cols;

	for (pt::u32 h_idx : ref_haps) {
		auto [col_idx, row_idx] = g.get_ref_gt_col_idx(h_idx);
		ns_cols.insert(col_idx);
		gt_cols[col_idx][row_idx] = "0";
	}

	// 1 because 0 is reserved for reference allele
	// pt::u32 i{1};
	for (auto alt_h_idx : alt_h_idxs) {
		auto [col_idx, row_idx] = g.get_ref_gt_col_idx(alt_h_idx);
		gt_cols[col_idx][row_idx] = std::to_string(1);
		ns_cols.insert(col_idx);
		// i++;
	}

	return {ns_cols.size(), gt_cols};
}

std::pair<pt::u32, std::vector<std::vector<std::string>>>
gen_gt_data(const bd::VG &g, const pga::minimal_rov &min_rov,
	    const pga::walk_to_alts_map &wta)
{
	std::vector<std::vector<std::string>> gt_cols =
		g.get_blank_genotype_cols();

	for (std::vector<std::string> &col : gt_cols)
		for (std::string &gt : col)
			gt = ".";

	// number of samples with data. cols that that are non-zero
	std::set<pt::u32> ns_cols;

	for (pt::u32 h_idx : min_rov.get_haps_matching_ref()) {
		auto [col_idx, row_idx] = g.get_ref_gt_col_idx(h_idx);
		ns_cols.insert(col_idx);
		gt_cols[col_idx][row_idx] = "0";
	}

	// 1 because 0 is reserved for reference allele
	pt::u32 i{1};
	for (const auto &[_, slices] : wta) {
		for (const pga::hap_slice &alt_as : slices) {
			auto [col_idx, row_idx] =
				g.get_ref_gt_col_idx(alt_as.ref_idx);
			gt_cols[col_idx][row_idx] = std::to_string(i);
			ns_cols.insert(col_idx);
		}
		i++;
	}

	return {ns_cols.size(), gt_cols};
}

void append_record(const bd::VG &g, pt::u32 ref_h_idx,
		   const pga::rov_boundaries &cxt,
		   const pvst::VertexBase *pvst_vtx_ptr, bool is_tangled,
		   const pga::minimal_rov &min_rov,
		   const pga::walk_to_alts_map &wta, pvr::var_type_e vt,
		   std::vector<VcfRec> &recs)
{
	if (wta.empty())
		return;

	pt::u32 pos = min_rov.get_ref_as().comp_pos(vt);

	std::set<pt::u32> ref_at_haps = min_rov.get_haps_matching_ref();

	VcfRec vcf_rec{ref_h_idx,
		       pos,
		       cxt.to_string(),
		       pvst_vtx_ptr->as_str(),
		       min_rov.get_ref_as(),
		       pvst_vtx_ptr->get_height(),
		       std::move(ref_at_haps),
		       vt,
		       is_tangled};

	for (const auto &[_, slices] : wta)
		vcf_rec.add_alt_set(slices);

	auto [ns, gt_data] = gen_gt_data(g, min_rov, wta);

	vcf_rec.add_gt_cols(std::move(gt_data));

	vcf_rec.set_ns(ns);

	recs.emplace_back(std::move(vcf_rec));
};

void gen_inv_recs(const bd::VG &g, const poi::it &it_,
		  std::vector<VcfRec> &recs)
{
	auto comp_alt_hap_slices =
		[&](const poi::vertex &v,
		    pt::u32 len) -> std::vector<pga::hap_slice>
	{
		const std::set<pt::u32> *alt_haps = v.get_len_alts(len);
		if (alt_haps == nullptr)
			return {};

		std::vector<pga::hap_slice> alt_set;
		for (pt::u32 alt_h_idx : *alt_haps)
			for (const poi::alt &a : v.get_alts(alt_h_idx))
				if (a.len == len)
					alt_set.emplace_back(
						g.get_ref_vec(a.h_idx)->walk,
						a.h_idx, a.h_start, len);

		return alt_set;
	};

	pt::u32 ref_h_idx = it_.get_ref_hap_idx();

	for (const auto &[_, v] : it_.get_vertices()) {
		pt::u32 ref_h_start = v.get_r_start();

		for (auto &[len, alts] : v.get_same_len_alts()) {

			const lq::ref_walk *ref_h_w =
				g.get_ref_vec(ref_h_idx)->walk;

			pt::u32 s = ref_h_start;
			pt::u32 t = ref_h_start + len - 1;

			std::string id =
				pv_cmp::format("|{}|{}|", ref_h_w->v_ids[s],
					       ref_h_w->v_ids[t]);

			pga::hap_slice ref_sl = {g.get_ref_vec(ref_h_idx)->walk,
						 ref_h_idx, ref_h_start, len};
			pt::u32 pos = ref_sl.comp_pos(pvr::var_type_e::inv);

			pt::u32 h = 0; // height

			std::set<pt::u32> ref_haps = {ref_h_idx};

			VcfRec vcf_rec{ref_h_idx,
				       pos,
				       id,
				       "",
				       ref_sl,
				       h,
				       std::move(ref_haps),
				       pvr::var_type_e::inv,
				       false};

			// TODO [A] actually fix this. Alts should always be
			// present for INV
			// if (alts.empty()) {
			//	WARN("1. No alt haplotypes for INV at pos {} "
			//	     " {}",
			//	     ref_sl.comp_pos(pvr::var_type_e::inv),
			//	     ref_h_start);
			//	continue;
			// }

			// std::vector<poi::alt> alts_ =
			//	v.get_alts(alt_h_idx, len);

			// for (const poi::alt &a : v.get_alts(alt_h_idx, len))
			// {
			// }

			// TODO [A] handle multiple alt haplotypes for
			// INV
			// std::vector<pga::hap_slice> alt_set;
			// for (auto alt_h_idx : alts) {
			//	for (const poi::alt &a :
			//	     v.get_alts(alt_h_idx)) {
			//		if (a.len == len) {
			//			alt_set.emplace_back(
			//				g.get_ref_vec(a.h_idx)
			//					->walk,
			//				a.h_idx, a.h_start,
			//				len);
			//		}
			//	}
			// }
			std::vector<pga::hap_slice> alt_set =
				comp_alt_hap_slices(v, len);

			if (alt_set.empty()) // no alts at that len
				continue;

			vcf_rec.add_alt_set(std::move(alt_set));

			auto [ns, gt_data] = gen_gt_data(g, {ref_h_idx}, alts);

			vcf_rec.add_gt_cols(std::move(gt_data));
			vcf_rec.set_ns(ns);

			recs.emplace_back(std::move(vcf_rec));
		}
	}
}

std::map<pt::idx_t, std::vector<VcfRec>> gen_exp_vcf_recs(const bd::VG &g,
							  const pga::trek &tk)
{

	std::map<pt::idx_t, std::vector<VcfRec>> tk_vcf_recs;

	const pvst::VertexBase *pvst_vtx_ptr = tk.get_pvst_vtx_const_ptr();

	for (pt::u32 ref_h_idx : tk.get_ref_haps()) {
		std::vector<VcfRec> recs;
		const pga::cxt_to_min_rov_map &d = tk.get_ref_recs(ref_h_idx);
		for (const auto &[cxt, min_rov] : d) {
			append_record(g, ref_h_idx, cxt, pvst_vtx_ptr,
				      tk.is_tangled(), min_rov,
				      min_rov.get_ins(), pvr::var_type_e::ins,
				      recs);

			append_record(g, ref_h_idx, cxt, pvst_vtx_ptr,
				      tk.is_tangled(), min_rov,
				      min_rov.get_subs(), pvr::var_type_e::sub,
				      recs);

			append_record(g, ref_h_idx, cxt, pvst_vtx_ptr,
				      tk.is_tangled(), min_rov,
				      min_rov.get_dels(), pvr::var_type_e::del,
				      recs);
		}

		tk_vcf_recs.insert({ref_h_idx, std::move(recs)});
	}

	return tk_vcf_recs;
}

void context_bound(const bd::VG &g, const std::vector<pga::trek> &treks,
		   VcfRecIdx &vcf_recs)
{
	for (pt::idx_t i{}; i < treks.size(); ++i) {
		const pga::trek &tk = treks[i];
		std::map<pt::idx_t, std::vector<VcfRec>> exp_vcf_recs =
			gen_exp_vcf_recs(g, tk);
		vcf_recs.ensure_append_recs(std::move(exp_vcf_recs));
	}
}

void context_free(const bd::VG &g, const std::vector<poi::it> &its,
		  VcfRecIdx &vcf_recs)
{
	for (const auto &i_tree : its) {
		pt::u32 ref_h_idx = i_tree.get_ref_hap_idx();
		auto &recs = vcf_recs.ensure_recs_mut(ref_h_idx);
		gen_inv_recs(g, i_tree, recs);
	}
}

VcfRecIdx gen_vcf_records(const bd::VG &g, const std::vector<pga::trek> &treks,
			  const std::vector<poi::it> &its)
{
	VcfRecIdx vcf_recs;

	context_bound(g, treks, vcf_recs);
	context_free(g, its, vcf_recs);

	// for (pt::idx_t i{}; i < treks.size(); ++i) {
	//	const pga::trek &tk = treks[i];
	//	std::map<pt::idx_t, std::vector<VcfRec>> exp_vcf_recs =
	//		gen_exp_vcf_recs(g, tk);
	//	vcf_recs.ensure_append_recs(std::move(exp_vcf_recs));
	// }

	return vcf_recs;
}

} // namespace povu::genomics::vcf
