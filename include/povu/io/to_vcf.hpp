#ifndef PV_IO_VCF_HPP
#define PV_IO_VCF_HPP

#include <fstream> // std::ofstream
#include <map>
#include <string>
#include <string_view>
#include <vector>

#include "common.hpp"
"#include "povu/common/app.hpp"
"#include "povu/common/compat.hpp"
"#include "povu/common/utils.hpp"
"#include "povu/genomics/vcf.hpp"
"#include "povu/graph/bidirected.hpp"

namespace povu::io::to_vcf
{
inline constexpr std::string_view MODULE = "povu::io::to_vcf";
namespace fs = std::filesystem;
namespace pgv = povu::genomics::vcf;
namespace bd = povu::bidirected;
namespace pu = povu::utils;

class VcfOutput
{
	// for combined output e.g to stdout
	std::ostream *combined_ = nullptr;

	// applies to split output
	// ref to ofs for the ref
	std::map<std::string, std::ofstream> files_;

	// keep default constructor private
	VcfOutput() = default;

public:
	// ---------------
	// factory methods
	// ---------------
	static VcfOutput to_stdout()
	{
		VcfOutput v;
		v.combined_ = &std::cout;
		return v;
	}

	static VcfOutput
	to_split_files(const std::vector<std::string> &ref_labels,
		       const fs::path &out_dir)
	{
		VcfOutput v;
		povu::io::common::create_dir_if_not_exists(out_dir);

		// open files for each ref label
		for (const auto &ref_label : ref_labels) {
			fs::path vcf_fp = out_dir / (ref_label + ".vcf");
			{ // if file exists, truncate
				std::ofstream(vcf_fp,
					      std::ios::out | std::ios::trunc);
			}
			std::ofstream ofs(vcf_fp,
					  std::ios::out | std::ios::app);
			if (!ofs) {
				ERR("Append open failed: {}",
				    fs::absolute(vcf_fp).string());
				std::exit(EXIT_FAILURE);
			}
			v.files_.emplace(ref_label, std::move(ofs));
		}

		return v;
	}

	// ----------------
	// getters
	// ----------------

	/**
	 * Get a stream for a given label (combined ignores label)
	 */
	std::ostream &stream_for(const std::string &ref_label)
	{
		if (combined_)
			return *combined_;

		if (pv_cmp::contains(files_, ref_label))
			return files_.at(ref_label);

		// TODO: [c] handle this before in the caller or at startup
		// try prefix match
		for (const auto &[k, _] : files_)
			if (pu::is_prefix(k, ref_label))
				return files_.at(k);

		throw std::runtime_error(
			"[VcfOutput::stream_for] Unknown label: " + ref_label);
	}

	/**
	 * Apply to every active stream (useful for common headers)
	 */
	void for_each_stream(const std::function<void(std::ostream &)> &fn)
	{
		if (combined_) {
			fn(*combined_);
			return;
		}

		for (auto &[_, ofs] : files_) {
			fn(ofs);
		}
		return;
	}

	void flush_all()
	{
		if (combined_) {
			combined_->flush();
			return;
		}

		for (auto &[_, ofs] : files_) {
			ofs.flush();
		}
		return;
	}
};

void init_vcfs(bd::VG &g, const std::vector<std::string> &sample_names,
	       VcfOutput &vout);
void write_vcfs(pgv::VcfRecIdx &vcf_recs, const bd::VG &g,
		const std::set<pt::id_t> &vcf_ref_ids, VcfOutput &vout,
		const core::config &app_config);
} // namespace povu::io::to_vcf

#endif // PV_IO_VCF_HPP
