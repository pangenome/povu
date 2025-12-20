#ifndef PV_IO_FROM_VCF_HPP
#define PV_IO_FROM_VCF_HPP

#include <filesystem> // for path
#include <map>	      // for map
#include <ostream>
#include <string>
#include <string_view> // for string_view
#include <vector>

#include "povu/common/constants.hpp" // for TAB_CHAR
#include "povu/common/core.hpp"	     // for pt
#include "povu/common/utils.hpp"     // for split

namespace povu::io::from_vcf
{
inline constexpr std::string_view MODULE = "povu::io::from_vcf";
namespace fs = std::filesystem;

struct at_meta {
	pt::u32 sample_idx;
	pt::u32 phase_idx;
};

struct gt_data {
private:
	std::map<pt::u32, std::vector<at_meta>> data; // allele_to_sample_phase

public:
	void insert(pt::u32 at_idx, const at_meta &ap)
	{
		this->data[at_idx].emplace_back(ap);
	}

	[[nodiscard]]
	const std::map<pt::u32, std::vector<at_meta>> &get_data() const
	{
		return this->data;
	}

	[[nodiscard]]
	std::vector<at_meta> get_meta_for_allele(pt::u32 at_idx) const
	{
		auto it = this->data.find(at_idx);
		if (it != this->data.end())
			return it->second;
		else
			return {};
	}
};

struct VCFRecord {
private:
	std::string chrom;
	pt::u32 pos;
	std::string id;
	std::string ref;
	std::vector<std::string> alts;

	std::vector<std::string> allele_traversals;
	gt_data genotypes;

	VCFRecord() = default;

	void extract_ats(const std::string &s)
	{
		std::vector<std::string> info_fields;
		povu::utils::split(s, ';', &info_fields);

		for (const auto &field : info_fields) {
			if (field.rfind("AT=", 0) == 0) {
				std::string ats_str = field.substr(3);
				povu::utils::split(ats_str, ',',
						   &this->allele_traversals);
			}
		}

		return;
	}

	void handle_gt(pt::u32 i, const std::string &gt_str, gt_data &a)
	{
		std::vector<std::string> d;
		povu::utils::split(gt_str, '|', &d);
		for (pt::u32 j = 0; j < d.size(); j++) {
			const std::string &allele = d[j];
			if (allele == ".")
				continue;

			// also allele index
			auto at_idx = static_cast<pt::u32>(std::stoul(allele));

			// i is the sample idx, j is the contig for that
			// sample or phase
			at_meta k{i, j};
			a.insert(at_idx, k);
		}
	}

	void handle_genotypes(const std::vector<std::string> &s,
			      pt::u32 start_col)
	{
		// gt_data a;

		// loop over all samples
		// genotype fields are from start_col to end
		for (pt::u32 i = start_col; i < s.size(); i++) {
			const std::string &gt_str = s[i];
			if (gt_str == ".")
				continue;

			this->handle_gt(i - start_col, gt_str, this->genotypes);
		}
	}

public:
	// --------------
	// factory method
	// --------------
	static VCFRecord from_row(const std::string &s)
	{
		VCFRecord rec;
		std::vector<std::string> fields;
		pu::split(s, pc::COL_SEP, &fields);

		rec.chrom = fields[0];
		rec.pos = static_cast<pt::u32>(std::stoul(fields[1]));
		rec.id = fields[2];
		rec.ref = fields[3];
		povu::utils::split(fields[4], ',', &rec.alts);

		rec.extract_ats(fields[7]);
		rec.handle_genotypes(fields, 9);
		return rec;
	}

	[[nodiscard]]
	const gt_data &get_genotypes() const
	{
		return this->genotypes;
	}

	[[nodiscard]]
	const std::string &get_at(pt::u32 at_idx) const
	{
		return this->allele_traversals[at_idx];
	}

	[[nodiscard]]
	pt::u32 get_pos() const
	{
		return this->pos;
	}

	[[nodiscard]]
	const std::string &get_id() const
	{
		return this->id;
	}

	void dbg_print(std::ostream &os) const
	{
		os << " chrom: " << this->chrom << "\n";
		os << " pos: " << this->pos << "\n";
		os << " id: " << this->id << "\n";
		os << " ref: " << this->ref << "\n";
		os << " alts: " << pu::concat_with(this->alts, ',') << "\n";

		auto N = static_cast<pt::u32>(this->allele_traversals.size());
		for (pt::u32 i = 0; i < N; i++) {
			os << allele_traversals[i] << "\t";
			std::vector<at_meta> d =
				this->genotypes.get_meta_for_allele(i);
			for (const auto &m : d) {
				os << "(" << m.sample_idx << "," << m.phase_idx
				   << ") ";
			}
			os << "\n";
		}
	}
};

struct VCFile {
private:
	pu::TwoWayMap<std::string, pt::u32> sns; // sample name to index
	std::vector<VCFRecord> records;

public:
	// -----------
	// constructor
	// -----------
	VCFile() = default;

	// -------
	// getters
	// -------

	[[nodiscard]]
	const std::vector<VCFRecord> &get_records() const
	{
		return this->records;
	}

	[[nodiscard]]
	pt::u32 get_sample_idx(const std::string &sample_name) const
	{
		return this->sns.get_value(sample_name);
	}

	[[nodiscard]]
	std::string get_sample_name(pt::u32 sample_idx) const
	{
		return this->sns.get_key(sample_idx);
	}

	// ---------
	// modifiers
	// ---------

	void add_sample_names(const std::string &s)
	{
		// split string by tabs
		std::vector<std::string> tokens;
		pu::split(s, '\t', &tokens);

		pt::u32 i{};
		for (; i < tokens.size(); i++)
			if (tokens[i] == "FORMAT")
				break;

		for (pt::u32 j = i + 1; j < tokens.size(); j++)
			this->sns.insert(tokens[j], j - (i + 1));
	}

	void dbg_print(std::ostream &os) const
	{
		pt::u32 i{};
		for (const auto &sn : this->sns.get_keys())
			os << i++ << " " << sn << ", ";

		os << "\n";

		for (const auto &rec : this->records)
			rec.dbg_print(os);
	}

	void add_record(VCFRecord &&rec)
	{
		this->records.emplace_back(rec);
	}
};

void read_vcf(const fs::path &fp, pt::u32 ll,
	      povu::io::from_vcf::VCFile &vcf_file);

} // namespace povu::io::from_vcf

#endif // PV_IO_FROM_VCF_HPP
