#include "zien/repeats/repeats.hpp"	  // for foo
#include "mto/from_vcf.hpp"		  // for VCFile
#include "povu/graph/bidirected.hpp"	  // for VG
#include "zien/components/components.hpp" // for Mode

namespace zien::components::repeats
{
/* repeats pane update function */
void update_repeats(const bd::VG &g, const mto::from_vcf::VCFile &vcf_file,
		    pt::u32 rec_idx, display_lines &pd)
{
	pd.reset();

	const mto::from_vcf::VCFRecord &rec =
		vcf_file.get_records().at(rec_idx);

	if (!rec.is_tangled())
		return;

	std::string ref_tag = rec.get_chrom();
	pt::u32 ref_h_id = *g.get_ref_id(ref_tag);
	pt::op_t<ptg::id_or_t> ef = rec.get_ef();
	auto [header_lens, hl_slices, lines] =
		zien::tui::repeats::foo(g, ref_h_id, ef, rec.get_pos());

	pt::u32 longest_header{};

	for (pt::u32 i{}; i < lines.size(); i++) {
		line_metadata &lm = pd.meta[i];
		lm.ref_name_pos = header_lens[i];
		lm.at_str_slices.push_back(hl_slices[i]);
		if (header_lens[i] > longest_header)
			longest_header = header_lens[i];
		if (i > 0 && i % 2 == 0)
			pd.group_lines.insert(i);
	}

	pd.lh = longest_header;
	pd.lines.swap(lines);
};
} // namespace zien::components::repeats
