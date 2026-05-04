#include "mto/from_vcf.hpp"

#include <set>
#include <string>
#include <vector>

#include "mto/common.hpp" // for fp_to_vector

#include <log.h>	   // for log_fatal
#include <quilt/shim.hpp>  // for contains
#include <quilt/types.hpp> // for qt
#include <quilt/utils.hpp> // for concat_with

namespace mto::from_vcf
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
			qt::u32 len_prefix = 17;
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
			log_fatal("Malformed VCF, Empty line found in VCF "
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
	if (!qs::contains(VCF_SUPPORTED_VERSIONS, vcf_version)) {
		std::string supported = pu::concat_with(
			std::vector<std::string>(VCF_SUPPORTED_VERSIONS.begin(),
						 VCF_SUPPORTED_VERSIONS.end()),
			',');

		log_fatal("Unsupported VCF version: %s. "
			  "Supported versions are: %s",
			  vcf_version.c_str(), supported.c_str());

		std::exit(EXIT_FAILURE);
	}

	return vcf_version;
}

pu::TwoWayMap<std::string, qt::u32> extract_sample_names(const std::string &s)
{
	// split string by tabs
	std::vector<std::string> tokens;
	pu::split(s, '\t', &tokens);

	qt::u32 i{};
	for (; i < tokens.size(); i++)
		if (tokens[i] == "FORMAT")
			break;

	pu::TwoWayMap<std::string, qt::u32> sn_to_idx;

	std::vector<std::string> sample_names;
	for (qt::u32 j = i + 1; j < tokens.size(); j++) {
		sample_names.push_back(tokens[j]);
		sn_to_idx.insert(tokens[j], j - (i + 1));
	}

	return sn_to_idx;
}

void read_vcf(const fs::path &fp, qt::u32 ll, mto::from_vcf::VCFile &vcf_file)
{
	if (ll > 1)
		log_info("Reading VCF from: {}", fp.string().c_str());

	std::vector<std::string> lines;
	mto::common::fp_to_vector(fp, &lines);

	std::vector<std::string> headers = extract_header_lines(lines);
	std::string vcf_version = check_version(headers);

	if (ll > 1)
		log_info("VCF version: {}", vcf_version.c_str());

	vcf_file.add_sample_names(headers.back());

	qt::u32 i{static_cast<qt::u32>(headers.size())};
	for (; i < lines.size(); i++) {
		const std::string &line = lines[i];
		auto rec = mto::from_vcf::VCFRecord::from_row(line);
		vcf_file.add_record(std::move(rec));
	}
}

// make this a utility function
std::vector<ptg::id_or_t> extract_v_id_ors(const std::string &at)
{
	std::vector<ptg::id_or_t> v_ids;

	ptg::or_e o;
	std::string v_id_str = "";

	auto update = [&v_id_str, &o, &v_ids]()
	{
		if (v_id_str != "") {
			auto v_id = static_cast<qt::id_t>(std::stoll(v_id_str));
			v_ids.emplace_back(ptg::id_or_t{v_id, o});
			v_id_str.clear();
		}
	};

	// zero is a > or <
	for (qt::u32 i{1}; i < at.size(); i++) {
		char c = at[i];

		if ((c == '>' || c == '<') && v_id_str != "")
			update();

		if (c == '>')
			o = ptg::or_e::forward;
		else if (c == '<')
			o = ptg::or_e::reverse;
		else if (std::isdigit(c))
			v_id_str += c;
	}

	update(); // final update

	return v_ids;
};

} // namespace mto::from_vcf
