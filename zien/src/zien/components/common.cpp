#include <liteseq/refs.h> // for ref_walk, ref

#include "mto/from_vcf.hpp"		  // for VCFile
#include "povu/graph/bidirected.hpp"	  // for VG
#include "zien/common/common.hpp"	  // for comp_update_refs
#include "zien/components/components.hpp" // for status_bar

#include <povu/refs/refs.hpp> // for ref_format_e, Ref, lq_strand_to_char

#include <quilt/types.hpp> // for qt

namespace zien::components::common
{
namespace lq = liteseq;
constexpr oza::refs::ref_format_e PN = oza::refs::ref_format_e::PANSN;
using namespace zien::components;

qt::id_t extract_anchor_v_id(const std::string &at)
{
	std::string v_id_str = "";
	// zero is a > or <
	for (qt::u32 i{1}; i < at.size(); i++) {
		char c = at[i];
		if (c == '>' || c == '<')
			break;

		if (std::isdigit(c))
			v_id_str += c;
	}

	return static_cast<qt::id_t>(std::stoll(v_id_str));
}

qt::u32 at_str_step_count(const std::string &at_str)
{
	qt::u32 count{};
	for (qt::u32 i{}; i < at_str.size(); i++) {
		char c = at_str[i];
		if (c == '>' || c == '<')
			count++;
	}

	return count;
}

std::set<qt::id_t> get_ref_ids(const bd::VG &g, const std::string &sn,
			       qt::u32 phase_idx)
{
	std::set<qt::id_t> filtered_ref_ids;

	for (qt::id_t r_id : g.get_refs_in_sample(sn)) {
		const oza::refs::Ref &r = g.get_ref_by_id(r_id);

		// TODO: find a better way to handle non PANSN
		if (r.get_format() != PN) // just trust it
			filtered_ref_ids.insert(r_id);
		else if (r.get_format() == PN && r.get_hap_id() == phase_idx)
			filtered_ref_ids.insert(r_id);
	}

	return filtered_ref_ids;
}

void comp_update_refs(const bd::VG &g, const mto::from_vcf::VCFile &vcf_file,
		      qt::u32 selected_rec_idx, qt::u32 at_idx,
		      display_lines &pd)
{
	qt::u32 line_count = pd.lines.size();
	qt::u32 ctr{line_count}; // counter for line metadata

	// what we compute and return
	// lm &meta = pd.meta;
	std::vector<std::string> &ref_lines = pd.lines;

	const mto::from_vcf::VCFRecord &rec =
		vcf_file.get_records().at(selected_rec_idx);
	const mto::from_vcf::gt_data &d = rec.get_genotypes();

	try {
		const std::vector<mto::from_vcf::at_meta> at_meta =
			d.get_data().at(at_idx);
	}
	catch (std::out_of_range &e) {
		std::cerr << " at idx " << at_idx << " rec idx "
			  << selected_rec_idx << "\n";
		std::exit(1);
	}

	const std::vector<mto::from_vcf::at_meta> at_meta =
		d.get_data().at(at_idx);

	const std::string &s = rec.get_at(at_idx);
	qt::u32 at_str_len = s.length();
	qt::u32 at_str_sc = at_str_step_count(s);
	qt::id_t anchor_v_id = extract_anchor_v_id(s);

	for (const auto &[sample_idx, phase_idx] : at_meta) {
		std::string sn = vcf_file.get_sample_name(sample_idx);
		std::string curr_l = "";

		std::set<qt::id_t> ref_ids = zien::common::get_ref_ids(
			g, vcf_file, sample_idx, phase_idx);

		for (qt::u32 h_idx : ref_ids) {

			const lq::ref_walk *rw = g.get_ref_vec(h_idx)->walk;
			qt::u32 N = rw->step_count;

			curr_l.clear();

			// curr_l += "[" + std::to_string(N) + "] ";
			curr_l += g.get_tag(h_idx);

			qt::u32 header_len = curr_l.length();

			if (header_len > pd.lh)
				pd.lh = header_len;

			// std::cerr << "Anchor v_id: " << anchor_v_id
			// << "\t"
			//	  << sn << "\n";

			qt::u32 anchor_v_idx = g.v_id_to_idx(anchor_v_id);
			const std::vector<qt::idx_t> &starts =
				g.get_vertex_ref_idxs(anchor_v_idx, h_idx);

			if (starts.empty())
				continue;

			line_metadata &lm = pd.meta[ctr]; // create if
							  // not exists
			lm.ref_name_pos = header_len;

			for (qt::u32 s{}; s < starts.size(); s++) {
				qt::u32 ref_pos = starts[s];
				qt::u32 i = ref_pos > 10 ? ref_pos - 10 : 0;
				qt::u32 end =
					std::min(ref_pos + at_str_sc + 10, N);

				for (; i < end; i++) {
					if (i == ref_pos) {
						qt::slice w{
							(qt::u32)curr_l.size(),
							at_str_len};

						lm.at_str_slices.emplace_back(
							w);
					}

					qt::id_t v_id = rw->v_ids[i];
					char o = oza::refs::lq_strand_to_char(
						rw->strands[i]);
					// char o = lq_strand_to_or_e(
					//	rw->strands[i]);

					curr_l += o;
					curr_l += std::to_string(v_id);
				}

				if (s + 1 < starts.size())
					curr_l += " ... ";
			}

			ref_lines.push_back(curr_l);
			ctr++;
			// pd.lines.push_back(curr_l);
			//  pd.p->lines.push_back(curr_l);
		}
	}

	return;
};

} // namespace zien::components::common
