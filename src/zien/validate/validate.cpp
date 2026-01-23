#include "zien/validate/validate.hpp"

#include <cctype>
#include <fstream>
#include <string>

#include <liteseq/refs.h> // for ref_walk, ref
#include <vector>

#include "mto/from_vcf.hpp" // for read_vcf

#include "povu/common/core.hpp"
#include "povu/graph/bidirected.hpp" // for bidirected
#include "povu/graph/types.hpp"	     // for or_e
#include "povu/refs/refs.hpp"	     // for Ref
#include "zien/common/common.hpp"    // for get_ref_ids

namespace zien::validate
{
namespace lq = liteseq;

constexpr povu::refs::ref_format_e PN = povu::refs::ref_format_e::PANSN;

ptg::walk_t at_to_walk(const std::string &at)
{
	ptg::walk_t w;
	std::string v_id_str = "";
	ptg::or_e o = ptg::or_e::forward;

	auto update_walk = [&]()
	{
		if (v_id_str.empty())
			return;

		ptg::id_or_t step;
		step.v_id = std::stoul(v_id_str);
		step.orientation = o;
		w.push_back(step);

		v_id_str = "";
		o = ptg::or_e::forward;
	};

	for (auto c : at) {
		if (c == '>' || c == '<') {
			update_walk();
			o = (c == '>') ? ptg::or_e::forward
				       : ptg::or_e::reverse;
			continue;
		}

		if (std::isdigit(c)) {
			v_id_str += c;
			continue;
		}
	}

	// last step
	update_walk();

	return w;
}

ptg::or_e lq_strand_to_or_e(lq::strand s)
{
	return (s == lq::strand::STRAND_FWD) ? ptg::or_e::forward
					     : ptg::or_e::reverse;
}

bool check_at(const bd::VG &g, const std::string &at, pt::idx_t ref_id)
{
	if (at.empty())
		return false;

	ptg::walk_t w = at_to_walk(at);
	pt::u32 at_len = w.size();
	pt::u32 s = g.v_id_to_idx(w.front().v_id);
	const std::vector<pt::idx_t> &positions =
		g.get_vertex_ref_idxs(s, ref_id);

	if (positions.empty())
		return false;

	const lq::ref_walk *rw = g.get_ref_vec(ref_id)->walk;

	bool any_of = false;

	auto check_at_loop = [at_len, rw, w](pt::u32 start_idx) -> bool
	{
		for (pt::u32 i{}; i < at_len; i++) {
			pt::u32 si = start_idx + i;
			pt::id_t ref_v_id = rw->v_ids[si];
			ptg::or_e ref_o = lq_strand_to_or_e(rw->strands[si]);

			ptg::id_or_t ref_step{ref_v_id, ref_o};

			if (ref_step != w[i])
				return false;
		}

		return true;
	};

	for (pt::u32 start_idx : positions)
		any_of = any_of || check_at_loop(start_idx);

	// we should not reach here
	return any_of;
}

std::set<pt::id_t> get_ref_ids_phased(const bd::VG &g, const std::string &sn,
				      pt::u32 phase_idx)
{
	// std::string sn = vcf_file.get_sample_name(sample_idx);
	std::set<pt::id_t> ref_ids = g.get_refs_in_sample(sn);

	std::set<pt::id_t> filtered_ref_ids;
	for (pt::id_t r_id : ref_ids) {
		const povu::refs::Ref &r = g.get_ref_by_id(r_id);

		// TODO: find a better way to handle non PANSN
		if (r.get_format() != PN) // just trust it
			filtered_ref_ids.insert(r_id);
		else if (r.get_format() == PN && r.get_hap_id() == phase_idx)
			filtered_ref_ids.insert(r_id);
	}

	return filtered_ref_ids;
}

bool any_contig(const bd::VG &g, const std::string &at,
		const std::set<pt::id_t> &ref_ids)
{
	for (pt::id_t ref_id : ref_ids)
		if (check_at(g, at, ref_id))
			return true;

	return false;
}

bool validate_rec(const bd::VG &g, const mto::from_vcf::VCFile &vcf_file,
		  pt::u32 rec_idx, std::ofstream *report_f, pt::u32 &err_recs)
{
	const mto::from_vcf::VCFRecord &rec = vcf_file.get_records()[rec_idx];

	// bool dbg = rec.get_pos() == 10318131 ? true : false;
	// if (dbg)
	//	rec.dbg_print(std::cerr);

	const mto::from_vcf::gt_data &d = rec.get_genotypes();
	for (const auto &[at_idx, at_meta] : d.get_data()) {
		const std::string &at = rec.get_at(at_idx);

		// std::cerr << at_idx << "\n";
		// std::cerr << at << "\n";

		for (const auto &[sample_idx, phase_idx] : at_meta) {

			std::set<pt::id_t> ref_ids = zien::common::get_ref_ids(
				g, vcf_file, sample_idx, phase_idx);

			std::string sn = vcf_file.get_sample_name(sample_idx);

			// // std::cerr << "sn " << sn
			// //	  << " Sample idx: " << sample_idx
			// //	  << " Phase idx: " << phase_idx << "\n";

			pt::u32 ploidy_id = g.get_ploidy_id(sn, phase_idx);

			// // ploidy is never 0, the else is always >1
			// std::set<pt::id_t> ref_ids =
			//	(g.get_ploidy(sn) == 0)
			//		? g.get_refs_in_sample(sn)
			//		: get_ref_ids_phased(g, sn, ploidy_id);

			// std::cerr << "Filtered Ref IDs: ["
			//	  << pu::concat_with(ref_ids, ',') << "]\n";

			// if (dbg) {
			//	std::cerr << "Filtered Ref IDs\n";
			//	std::cerr << pu::concat_with(ref_ids, ',')
			//		  << "\n";
			// }

			if (!any_contig(g, at, ref_ids)) {
				// std::string sn =
				//	vcf_file.get_sample_name(sample_idx);
				// ERR("Validation failed.\n "
				//     "Record {} id {} pos {} sample {} AT {}
				//     {} "
				//     "\n",
				//     rec_idx, rec.get_id(), rec.get_pos(), sn,
				//     at, h_id);
				if (report_f != nullptr)
					*report_f << pv_cmp::format(
						"{}\t{}\t{}\t{}\t{}\t{}\n",
						rec_idx, rec.get_id(),
						rec.get_pos(), at, ploidy_id,
						sn);

				err_recs++;
				// std::cout << rec_idx << "\t" << rec.get_id()
				//	  << "\t" << rec.get_pos() << "\t" << at
				//	  << "\t" << h_id << "\t" << sn << "\n";
				return false;
			}
		}
	}

	return true;
}

void write_summary(const core::config &app_config, pt::u32 err_recs, pt::u32 N)
{
	std::string summary_fp = pv_cmp::format(
		"{}/summary.txt",
		std::string{app_config.get_output_dir()}); // file path and name

	std::ofstream os(summary_fp);
	if (!os.is_open()) {
		PL_ERR("Could not open file {}", summary_fp);
		std::exit(EXIT_FAILURE);
	}

	if (err_recs == 0) {
		os << "No errors\n";
		return;
	}

	// comp err records as a percentage
	double fail_rate = (static_cast<double>(err_recs) / N) * 100;
	// int rounded_fail_rate = static_cast<int>(
	//	std::round(fail_rate)); // Round to nearest integer

	os << "Error rate: " << std::fixed << std::setprecision(2) << fail_rate
	   << "%\n";
}

std::vector<pt::u32> validate_vcf_records(const bd::VG &g,
					  const mto::from_vcf::VCFile &vcf_file,
					  const core::config &app_config,
					  bool output_to_file)
{
	std::ofstream file_stream; // Lives on the stack
	std::ofstream *report_os = nullptr;

	if (output_to_file) {
		// file path & name
		std::string report_fp = pv_cmp::format(
			"{}/report.tsv",
			std::string{app_config.get_output_dir()});

		file_stream.open(report_fp);
		if (!file_stream.is_open()) {
			PL_ERR("Could not open file {}", report_fp);
			std::exit(EXIT_FAILURE);
		}

		report_os = &file_stream;
		*report_os << "rec_idx\tID\tPOS\tAT\tHap ID\tsample\n";
	}

	std::vector<pt::u32> invalid_rec_indices;

	pt::u32 err_recs{};
	const pt::u32 N = vcf_file.get_records().size();
	for (pt::u32 rec_idx{}; rec_idx < N; rec_idx++)
		if (!validate_rec(g, vcf_file, rec_idx, report_os, err_recs))
			invalid_rec_indices.emplace_back(rec_idx);

	write_summary(app_config, err_recs, N);

	return invalid_rec_indices;
}
} // namespace zien::validate
