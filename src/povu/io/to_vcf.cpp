#include "povu/io/to_vcf.hpp"

#include <sstream> // for basic_ostringstream
#include <sys/types.h>

#include "fmt/core.h" // for format

#include "ita/variation/rov.hpp" // for var_type_e
				 //
#include "povu/common/core.hpp"	 // for pt
#include "povu/refs/refs.hpp"	 // for Ref, pr

namespace povu::io::to_vcf
{
namespace pgt = povu::types::graph;

constexpr std::string_view VCF_VERSION = "4.2";

void write_header_common(std::ostream &os)
{
	// clang-format off
	os << pv_cmp::format("##fileformat=VCFv{}\n", VCF_VERSION);
	os << pv_cmp::format("##fileDate={}\n", pu::today());
	os << "##source=povu\n";
	os << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	os << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">\n";
	os << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele traversal path through the graph\">\n";
	os << "##INFO=<ID=AN,Number=1,Type=String,Description=\"Total number of alleles in called genotypes\">\n";
	os << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency in the population\">\n";
	os << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">\n";
	os << "##INFO=<ID=VARTYPE,Number=1,Type=String,Description=\"Type of variation: INS (insertion), DEL (deletion), SUB (substitution), INV(inversion) \">\n";
	os << "##INFO=<ID=TANGLED,Number=1,Type=String,Description=\"Variant lies in a tangled region of the graph: T or F\">\n";
	os << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the PVST (0=top level)\">\n";
	os << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	// clang-format on
}

void write_header_contig_line(const pr::Ref &r, std::ostream &os)
{
	os << pv_cmp::format("##contig=<ID={},length={}>\n", r.tag(),
			     r.get_length());
}

void write_col_header(const std::vector<std::string> &genotype_col_names,
		      std::ostream &os)
{
	os << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT";
	const char sep = '\t';
	for (const auto &col_name : genotype_col_names)
		os << sep << col_name;
	os << "\n";
}

void init_vcfs(bd::VG &g, const std::vector<std::string> &ref_name_prefixes,
	       VcfOutput &vout)
{
	// write common header lines
	vout.for_each_stream([&](std::ostream &os)
			     { write_header_common(os); });

	// add contig lines
	for (const auto &rn_pref : ref_name_prefixes) {
		std::ostream &os = vout.stream_for_ref_label(rn_pref);
		std::set<pt::id_t> ref_ids = g.get_refs_in_sample(rn_pref);
		for (pt::id_t ref_id : ref_ids) {
			const pr::Ref &ref = g.get_ref_by_id(ref_id);
			write_header_contig_line(ref, os);
		}
	}

	// write column header line
	vout.for_each_stream(
		[&](std::ostream &os)
		{ write_col_header(g.get_genotype_col_names(), os); });

	vout.flush_all();

	return;
}

/**
 * @brief write a VCF record to output stream
 */
void write_rec(const bd::VG &g, iv::VcfRec &r, const std::string &chrom,
	       std::ostream &os)
{
	std::ostringstream info_field;
	// clang-format off
	info_field << "AC=" << r.get_ac()
		   << ";AF=" << r.get_af()
		   << ";AN=" << r.get_an()
		   << ";NS=" << r.get_ns()
		   << ";AT=" << r.get_at()
		   << ";VARTYPE=" << ir::to_string_view(r.get_var_type())
		   << ";TANGLED=" << (r.is_tangled() ? "T" : "F");
	// clang-format on

	if (r.get_var_type() != ir::var_type_e::inv) {
		// clang-format off
		info_field << ";ES=" << r.get_enc_flubble()
			   << ";LV=" << (r.get_height() - 1);
		// clang-format on
	}

	// clang-format off
	os << chrom << "\t"
	   << r.get_pos() << "\t"
	   << r.get_id() << "\t"
	   << r.get_ref_as_dna_str(g) << "\t"
	   << r.get_alts_as_dna_str(g) << "\t"
	   << r.get_qual() << "\t"
	   << r.get_filter() << "\t"
	   << info_field.str() << "\t"
	   << r.get_format() << "\t"
	   << r.get_genotype_fields() << "\n";
	// clang-format on

	return;
}

void write_vcfs(iv::VcfRecIdx &vcf_recs, const bd::VG &g, VcfOutput &vout,
		const core::config &app_config)
{
	// Cache stdout stream once to avoid repeatedly asking for it.
	const bool to_stdout = app_config.get_stdout_vcf();
	std::ostream *stdout_os =
		to_stdout ? &vout.stream_for_combined() : nullptr;

	for (auto &[ref_id, recs] : vcf_recs.get_recs_mut()) {
		std::string ref_tag = g.get_ref_by_id(ref_id).tag();
		std::ostream &os =
			to_stdout ? *stdout_os : vout.stream_for_ref_id(ref_id);

		for (iv::VcfRec &r : recs)
			write_rec(g, r, ref_tag, os);
	}

	return;
}

} // namespace povu::io::to_vcf
