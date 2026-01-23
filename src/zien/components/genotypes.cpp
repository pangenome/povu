#include <ncurses.h>

#include "mto/from_vcf.hpp" // for VCFile
#include "povu/common/core.hpp"
#include "povu/graph/bidirected.hpp"	  // for VG
#include "zien/components/components.hpp" // for Mode

namespace zien::components::genotypes
{

std::vector<std::string> comp_gt_data(const bd::VG &g,
				      const mto::from_vcf::VCFile &vcf_file,
				      pt::u32 selected_rec_idx)
{
	const mto::from_vcf::VCFRecord &rec =
		vcf_file.get_records().at(selected_rec_idx);
	const mto::from_vcf::gt_data &d = rec.get_genotypes();

	pt::u32 at_count = rec.get_at_count();

	std::vector<std::string> hap_lines;
	std::string hl;

	for (pt::u32 at_idx{}; at_idx < at_count; at_idx++) {
		const std::string &at = rec.get_at(at_idx);
		hl += at;
		hl += "\t";
	}

	hap_lines.emplace_back(hl);
	hl.clear();

	// count rows
	pt::u32 row_count{};
	for (pt::u32 at_idx{}; at_idx < at_count; at_idx++) {
		if (!pv_cmp::contains(d.get_data(), at_idx)) {
			continue;
		}
		// try {
		//	const std::vector<povu::io::from_vcf::at_meta>
		// at_meta =		d.get_data().at(at_idx);
		// }
		// catch (const std::out_of_range &e) {
		//	continue;
		// }

		const std::vector<mto::from_vcf::at_meta> &at_meta =
			d.get_data().at(at_idx);
		pt::u32 N = at_meta.size();
		if (N > row_count)
			row_count = N;
	}

	for (pt::u32 row_idx{}; row_idx < row_count; row_idx++) {
		for (pt::u32 at_idx{}; at_idx < at_count; at_idx++) {

			if (!pv_cmp::contains(d.get_data(), at_idx)) {
				// AT is not supported by any haps
				hl += "?";
				hl += "\t";
				continue;
			}

			// try {
			//	const
			// std::vector<povu::io::from_vcf::at_meta>
			//		at_meta =
			// d.get_data().at(at_idx);
			// }
			// catch (const std::out_of_range &e) {
			//	hl += "\t";
			//	continue;
			// }

			const std::vector<mto::from_vcf::at_meta> &at_meta =
				d.get_data().at(at_idx);

			if (row_idx < at_meta.size()) {
				const auto &[sample_idx, phase_idx] =
					at_meta[row_idx];
				std::string sn =
					vcf_file.get_sample_name(sample_idx);

				std::string name; // TODO: rename to row data or
						  // row string

				if (g.get_ploidy(sn) == 0) {
					pt::u32 hap_id =
						*g.get_refs_in_sample(sn)
							 .begin();
					name = g.get_tag(hap_id);
				}
				else {
					pt::u32 ploidy_id =
						g.get_ploidy_id(sn, phase_idx);

					name = sn + '#' +
					       std::to_string(ploidy_id);
				}

				hl += name;
				hl += "\t";
			}
			else {
				hl += "\t";
			}
		}
		hap_lines.emplace_back(hl);
		hl.clear();
	}

	return hap_lines;
};

/* top-right (haplotype) pane */
void update_haps(const bd::VG &g, const mto::from_vcf::VCFile &vcf_file,
		 pt::u32 selected_rec_idx, display_lines &pd)
{
	std::vector<std::string> hap_lines =
		comp_gt_data(g, vcf_file, selected_rec_idx);

	for (const std::string &hl : hap_lines)
		pd.lines.push_back(hl);
};

} // namespace zien::components::genotypes
