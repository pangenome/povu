#include "povu/io/from_vcf.hpp"

// #include <iostream> // for std::cerr, std::cout
// #include <map>	    // for map
#include <set>
#include <string>
#include <vector>

#include "povu/common/compat.hpp" // for pv_cmp
#include "povu/common/core.hpp"	  // for pt
#include "povu/common/log.hpp"	  // for INFO
#include "povu/common/utils.hpp"  // for concat_with
#include "povu/io/common.hpp"	  // for fp_to_vector

namespace povu::io::from_vcf
{

const std::set<std::string> VCF_SUPPORTED_VERSIONS = {
	"4.2",
};

std::string get_version(const std::vector<std::string> &header_lines)
{
	std::string vcf_version = "unknown";

	for (const auto &line : header_lines) {
		if (line.rfind("##fileformat=", 0) == 0) {
			// length of "##fileformat=VCFv" is 17
			pt::u32 len_prefix = 17;
			std::string version = line.substr(
				len_prefix, line.size() - len_prefix);
			vcf_version = version;
			break;
		}
	}

	return vcf_version;
}

std::vector<std::string>
extract_header_lines(const std::vector<std::string> &all_lines)
{
	std::vector<std::string> header_lines;

	for (const std::string &line : all_lines) {
		if (line.empty()) {
			ERR("Malformed VCF, Empty line found in VCF "
			    "header.");
			std::exit(EXIT_FAILURE);
		}

		if (line[0] == '#')
			header_lines.push_back(line);
		else
			return header_lines;
	}

	return header_lines;
}

std::string check_version(const std::vector<std::string> &headers)
{
	std::string vcf_version = get_version(headers);
	if (!pv_cmp::contains(VCF_SUPPORTED_VERSIONS, vcf_version)) {
		std::string supported = pu::concat_with(
			std::vector<std::string>(VCF_SUPPORTED_VERSIONS.begin(),
						 VCF_SUPPORTED_VERSIONS.end()),
			',');

		ERR("Unsupported VCF version: {}. Supported versions are: {}",
		    vcf_version, supported);

		std::exit(EXIT_FAILURE);
	}

	return vcf_version;
}

pu::TwoWayMap<std::string, pt::u32> extract_sample_names(const std::string &s)
{
	// split string by tabs
	std::vector<std::string> tokens;
	pu::split(s, '\t', &tokens);

	pt::u32 i{};
	for (; i < tokens.size(); i++)
		if (tokens[i] == "FORMAT")
			break;

	pu::TwoWayMap<std::string, pt::u32> sn_to_idx;

	std::vector<std::string> sample_names;
	for (pt::u32 j = i + 1; j < tokens.size(); j++) {
		sample_names.push_back(tokens[j]);
		sn_to_idx.insert(tokens[j], j - (i + 1));
	}

	return sn_to_idx;
}

void read_vcf(const fs::path &fp, pt::u32 ll,
	      povu::io::from_vcf::VCFile &vcf_file)
{
	if (ll > 1)
		INFO("Reading VCF from: {}", fp.string());

	std::vector<std::string> lines;
	povu::io::common::fp_to_vector(fp, &lines);

	std::vector<std::string> headers = extract_header_lines(lines);
	std::string vcf_version = check_version(headers);

	if (ll > 1)
		INFO("VCF version: {}", vcf_version);

	// pu::TwoWayMap<std::string, pt::u32> sns =
	//	extract_sample_names(headers.back());
	// INFO("Sample names ({}): {}", sns.size(), pu::concat_with(sns, ','));

	// povu::io::from_vcf::VCFile vcf_file;
	vcf_file.add_sample_names(headers.back());

	pt::u32 i{static_cast<pt::u32>(headers.size())};
	for (; i < lines.size(); i++) {
		const std::string &line = lines[i];
		auto rec = povu::io::from_vcf::VCFRecord::from_row(line);
		vcf_file.add_record(std::move(rec));
		// break;
	}

	// vcf_file.dbg_print(std::cerr);
	// return;
}

} // namespace povu::io::from_vcf
