#include "povu/io/to_vcf.hpp"

#include <sstream> // for basic_ostringstream
#include <sys/types.h>

#include "fmt/core.h" // for format
#include "povu/common/core.hpp"
#include "povu/genomics/allele.hpp" // for allele_slice_t
#include "povu/graph/types.hpp"	    // for step_t, or_e
#include "povu/refs/refs.hpp"	    // for Ref, pr
#include "povu/variation/rov.hpp"   // for var_type_e

namespace povu::io::to_vcf
{
namespace pgt = povu::types::graph;
namespace pga = povu::genomics::allele;

// clang-format off
void write_header_common(std::ostream &os) {
	os << "##fileformat=VCFv4.2\n";
	os << pv_cmp::format("##fileDate={}\n", pu::today());
	os << "##source=povu\n";
	os << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
	os << "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Total number of alternate alleles in called genotypes\">\n";
	os << "##INFO=<ID=AT,Number=R,Type=String,Description=\"Allele traversal path through the graph\">\n";
	os << "##INFO=<ID=AN,Number=1,Type=String,Description=\"Total number of alleles in called genotypes\">\n";
	os << "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency in the population\">\n";
	os << "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of samples with data\">\n";
	os << "##INFO=<ID=VARTYPE,Number=1,Type=String,Description=\"Type of variation: INS (insertion), DEL (deletion), SUB (substitution)\">\n";
	os << "##INFO=<ID=TANGLED,Number=1,Type=String,Description=\"Variant lies in a tangled region of the graph: T or F\">\n";
	os << "##INFO=<ID=LV,Number=1,Type=Integer,Description=\"Level in the PVST (0=top level)\">\n";
	os << "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n";
}

// clang-format on

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
void write_rec(const bd::VG &g, pgv::VcfRec &r, const std::string &chrom,
	       std::ostream &os)
{
	// std::string s = ">185>188";
	// bool dbg = (r.get_id() == s) ? true : false;
	// dbg = false;

	pvr::var_type_e var_typ = r.get_var_type();

	r.gen_rec_data_lookups(g); // ensure lookups are generated
	const pt::idx_t REF_AT_IDX = 0;
	std::vector<pt::idx_t> alts = r.get_unique_alt_idxs();

	auto get_label = [&](const pgt::step_t &s) -> std::string
	{
		auto [v_id, o] = s;
		return (o == pgt::or_e::forward)
			       ? g.get_vertex_by_id(v_id).get_label()
			       : g.get_vertex_by_id(v_id).get_rc_label();
	};

	auto slice_to_dna_str =
		[&](const pga::allele_slice_t &as) -> std::string
	{
		std::string dna_str = "";
		bool is_fwd = as.get_or() == pgt::or_e::forward;

		pt::u32 i = as.ref_start_idx;
		pt::u32 N = is_fwd ? as.ref_start_idx + as.len
				   : as.ref_start_idx - as.len;

		// 1) Anchor base for deletions & insertions
		if (var_typ == pvr::var_type_e::del ||
		    var_typ == pvr::var_type_e::ins) {
			// grab the first step’s label, take its last character
			const pgt::step_t &s = as.get_step(i);
			auto lbl = get_label(s);
			dna_str.push_back(is_fwd ? lbl.back() : lbl.front());
		}

		// 2) Middle steps (for all types) does nothing for deletions

		int8_t step = is_fwd ? 1 : -1;
		pt::u32 start = is_fwd ? ++i : --i;
		pt::u32 end = is_fwd ? --N : ++N;

		for (pt::u32 j{start}; j != end; j += step) {
			const pgt::step_t &s = as.get_step(j);
			dna_str += get_label(s);
		}

		// 3) If nothing got added, use “.”
		return dna_str.empty() ? std::string{"."} : dna_str;
	};

	auto slices_to_dna_str =
		[&](const std::vector<pt::idx_t> &idxs) -> std::string
	{
		std::string dna_strs = "";
		std::string sep = "";
		for (auto idx : idxs) {
			const pga::allele_slice_t &as = r.get_slice(idx);
			dna_strs += sep;
			// convert to DNA string and takes care of indels
			dna_strs += slice_to_dna_str(as);
			sep = ",";
		}

		return dna_strs;
	};

	auto slices_as_at_str = [&](const std::vector<pt::idx_t> &idxs,
				    pvr::var_type_e vt) -> std::string
	{
		std::string str = "";
		std::string sep = ""; // no leading comma
		for (auto i : idxs) {
			str += sep;
			str += r.get_slice(i).as_str(var_typ);
			sep = ",";
		}
		return str;
	};

	auto allele_traversals = [&]() -> std::string
	{
		std::string s;
		s += slices_as_at_str(std::vector<pt::idx_t>{REF_AT_IDX},
				      var_typ);
		s += ",";
		s += slices_as_at_str(alts, pvr::covariant(var_typ));
		return s;
	};

	// if (dbg) {
	//	std::cerr << slices_as_at_str(
	//			     std::vector<pt::idx_t>{REF_AT_IDX},
	//			     var_typ)
	//		  << "\n";
	//	std::cerr << slices_as_at_str(alts, pvr::covariant(var_typ))
	//		  << "\n";
	// }

	// Both ref_dna and alt_dna are plain std::string values over the DNA
	// letters {A, C, G, T}.
	//   ref_dna  is a single contiguous sequence.
	//   alt_dna  may contain one or more such sequences, separated by
	//   commas.
	std::string ref_dna = slice_to_dna_str(r.get_slice(REF_AT_IDX));
	std::string alt_dna = slices_to_dna_str(alts);

	std::ostringstream info_field;
	// clang-format off
	info_field << "AC=" << r.get_ac()
		   << ";AF=" << r.get_af()
		   << ";AN=" << r.get_an()
		   << ";NS=" << r.get_ns()
		   << ";AT=" << allele_traversals()
		   << ";VARTYPE=" << pvr::to_string_view(var_typ)
		   << ";TANGLED=" << (r.is_tangled() ? "T" : "F")
		   << ";LV=" << (r.get_height() - 1);
	// clang-format on

	// clang-format off
	os << chrom << "\t"
	   << r.get_pos() << "\t"
	   << r.get_id() << "\t"
	   << ref_dna << "\t"
	   << alt_dna << "\t"
	   << r.get_qual() << "\t"
	   << r.get_filter() << "\t"
	   << info_field.str() << "\t"
	   << r.get_format() << "\t"
	   << r.get_genotype_fields() << "\n";
	// clang-format on

	return;
}

void write_vcfs(pgv::VcfRecIdx &vcf_recs, const bd::VG &g, VcfOutput &vout,
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

		for (pgv::VcfRec &r : recs)
			write_rec(g, r, ref_tag, os);
	}

	return;
}

} // namespace povu::io::to_vcf
