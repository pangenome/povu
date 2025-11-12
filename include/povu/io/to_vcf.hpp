#ifndef PV_IO_VCF_HPP
#define PV_IO_VCF_HPP

#include <cstddef>
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
	std::vector<std::ofstream> all_ofs_;
	std::map<std::string, pt::idx_t> label_to_ofs_idx_;
	std::map<pt::idx_t, pt::idx_t> ref_id_to_ofs_idx_;

	// keep default constructor private
	VcfOutput() = default;

	// ----------------
	// private methods
	// ----------------
	std::size_t create_ofs(const fs::path &fp)
	{
		{ // if file exists, truncate
			std::ofstream(fp, std::ios::out | std::ios::trunc);
		}
		std::ofstream ofs(fp, std::ios::out | std::ios::app);
		if (!ofs) {
			ERR("Append open failed: {}",
			    fs::absolute(fp).string());
			std::exit(EXIT_FAILURE);
		}
		this->all_ofs_.push_back(std::move(ofs));

		return this->all_ofs_.size() - 1;
	}

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

	/**
	 * s_to_r sample to ref_ids
	 */
	static VcfOutput
	to_split_files(const fs::path &out_dir,
		       const std::map<std::string, std::set<pt::id_t>> &s_to_r)
	{
		VcfOutput v;
		povu::io::common::create_dir_if_not_exists(out_dir);

		// open files for each ref label
		for (const auto &[bn, ref_ids] : s_to_r) {
			fs::path vcf_fp = out_dir / (bn + ".vcf");
			std::size_t ofs_idx = v.create_ofs(vcf_fp);
			v.label_to_ofs_idx_[bn] = ofs_idx;
			for (pt::id_t ref_id : ref_ids)
				v.ref_id_to_ofs_idx_[ref_id] = ofs_idx;
		}

		return v;
	}

	// -------
	// getters
	// -------

	// TODO: [c] merge stream_for queries?

	std::ostream &stream_for_combined()
	{
		if (!combined_)
			throw std::runtime_error(
				"[VcfOutput::stream_for_combined] Not a "
				"combined output");

		return *combined_;
	}

	/**
	 * Get a stream for a given label (combined ignores label)
	 */
	std::ostream &stream_for_ref_label(const std::string &ref_label)
	{
		if (combined_)
			return *combined_;

		if (pv_cmp::contains(label_to_ofs_idx_, ref_label))
			return this->all_ofs_[label_to_ofs_idx_[ref_label]];

		// TODO: [c] handle this before in the caller or at startup
		// try prefix match
		for (const auto &[k, _] : label_to_ofs_idx_)
			if (pu::is_prefix(k, ref_label))
				return this->all_ofs_[label_to_ofs_idx_[k]];

		throw std::runtime_error(
			"[VcfOutput::stream_for] Unknown label: " + ref_label);
	}

	std::ostream &stream_for_ref_id(pt::idx_t ref_id)
	{
		if (combined_)
			return *combined_;

		if (pv_cmp::contains(ref_id_to_ofs_idx_, ref_id))
			return this->all_ofs_[ref_id_to_ofs_idx_[ref_id]];

		throw std::runtime_error(
			"[VcfOutput::stream_for] Unknown ref id: " + ref_id);
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

		for (auto &ofs : all_ofs_)
			fn(ofs);

		return;
	}

	void flush_all()
	{
		if (combined_) {
			combined_->flush();
			return;
		}

		for (auto &ofs : all_ofs_)
			ofs.flush();

		return;
	}
};

void init_vcfs(bd::VG &g, const std::vector<std::string> &sample_names,
	       VcfOutput &vout);
void write_vcfs(pgv::VcfRecIdx &vcf_recs, const bd::VG &g, VcfOutput &vout,
		const core::config &app_config);
} // namespace povu::io::to_vcf

#endif // PV_IO_VCF_HPP
