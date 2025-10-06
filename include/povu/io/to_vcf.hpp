#ifndef PV_IO_VCF_HPP
#define PV_IO_VCF_HPP

#include <cstdlib>     // for exit, EXIT_FAILURE
#include <filesystem>  // for path, absolute, operator/
#include <fstream>     // for basic_ofstream, basic_ios, ios
#include <functional>  // for function
#include <iostream>    // for cout
#include <map>	       // for map, operator!=
#include <set>	       // for set
#include <stdexcept>   // for runtime_error
#include <string>      // for basic_string, operator+, operator<
#include <string_view> // for string_view
#include <utility>     // for get, move
#include <vector>      // for vector

#include "common.hpp"		     // for create_dir_if_not_exists
#include "povu/common/app.hpp"	     // for config
#include "povu/common/compat.hpp"    // for contains, pv_cmp
#include "povu/common/core.hpp"	     // for id_t, pt
#include "povu/common/log.hpp"	     // for ERR
#include "povu/common/utils.hpp"     // for is_prefix
#include "povu/genomics/vcf.hpp"     // for VcfRecIdx
#include "povu/graph/bidirected.hpp" // for VG

namespace povu::io::to_vcf
{
inline constexpr std::string_view MODULE = "povu::io::to_vcf";
namespace fs = std::filesystem;
namespace pgv = povu::genomics::vcf;
namespace bd = povu::bidirected;

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
