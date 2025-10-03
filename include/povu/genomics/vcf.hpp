#ifndef POVU_GENOMICS_VCF_HPP
#define POVU_GENOMICS_VCF_HPP

#include <algorithm>
#include <map>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

"#include "povu/common/compat.hpp"
"#include "povu/graph/bidirected.hpp"
"#include "povu/graph/types.hpp"
#include "allele.hpp"

namespace povu::genomics::vcf
{
inline constexpr std::string_view MODULE = "povu::genomics::vcf";

namespace bd = povu::bidirected;
namespace pgt = povu::types::graph;
namespace pvst = povu::pvst;
namespace pga = povu::genomics::allele;
namespace pgg = povu::genomics::graph;

enum class var_type_e {
	del, // deletion
	ins, // insertion
	sub, // substitution
	und  // undetermined
};

constexpr std::string_view to_string_view(var_type_e vt) noexcept
{
	switch (vt) {
	case var_type_e::del:
		return "DEL";
	case var_type_e::ins:
		return "INS";
	case var_type_e::sub:
		return "SUB";
	case var_type_e::und:
		return "UND";
	}
	// optional: handle out-of-range
	return "??";
}

inline std::ostream &operator<<(std::ostream &os, var_type_e vt)
{
	return os << to_string_view(vt);
}

class VcfRec
{
	pt::id_t ref_id_; // chrom TODO: does this still apply with tags?
	pt::idx_t pos_;	  // 1-based step idx
	std::string id_;  // start and end of a RoV e.g >1>4

	// allele slices
	// allele slice at idx 0 is always the ref allele
	// subsequent indices are alt alleles
	std::vector<pga::allele_slice_t> ats_;

	inline static const std::string qual = "60";	 // fixed at 60
	inline static const std::string filter = "PASS"; // fixed as pass
	inline static const std::string format = "GT";	 // fixed as pass

	pt::idx_t height_;    // height of the pvst node in the tree
	var_type_e var_type_; // type of the variant, e.g. del, ins, sub, und
	bool is_tangled_ = false; // is true when tangling exists, i.e. when a
				  // walk traverses an RoV more than once

	/* info */
	// number of refs in a given walk
	std::vector<pt::idx_t> ref_count;
	std::vector<pt::idx_t> allele_counts_; // count for each allele (ref at
					       // idx 0, alts at idx 1+)
	std::vector<double> af_; // allele frequency for each alt allele
	pt::idx_t an_ = 0;	 // total number of alleles in called genotypes
	// pt::idx_t ns_ = 0;       // number of samples with data
	std::map<pt::idx_t, pt::idx_t>
		col_to_allele_count_; // genotype column idx to allele count

	// genotype
	std::vector<std::vector<std::string>> genotype_cols_;

	// lookups
	// each value represents a unique alt allele idx
	std::vector<pt::idx_t> unique_alts_;
	// maps alt allele idx (0-based) to unique alt allele idx (1-based)
	std::vector<pt::idx_t> alt_at_to_unique_alt_idx_;

	// -----------------
	// private method(s)
	// ------------------

	void comp_unique_alt_ats()
	{
		std::map<std::tuple<pt::idx_t, pt::idx_t, pt::idx_t>, pt::idx_t>
			seen;
		for (pt::idx_t i{}; i < this->ats_.size(); ++i) {
			const auto &alt_at = this->ats_.at(i);
			auto key = std::make_tuple(alt_at.walk_idx,
						   alt_at.walk_start_idx,
						   alt_at.len);
			if (!pv_cmp::contains(seen, key)) {
				this->unique_alts_.emplace_back(i);
				seen[key] = this->unique_alts_.size() - 1;
			}

			this->alt_at_to_unique_alt_idx_.emplace_back(seen[key]);
		}
		// Initialize allele_counts_ with the right size (ref + unique
		// alts)
		this->allele_counts_.resize(this->unique_alts_.size(), 0);
	}

	void set_genotype_data(const bd::VG &g)
	{
		// Initialize all genotypes to 0 (reference)
		for (auto &col : this->genotype_cols_) {
			for (auto &gt : col) {
				gt = "0";
			}
		}

		for (pt::idx_t i{}; i < this->ats_.size(); ++i) {
			const pga::allele_slice_t &sl = this->ats_.at(i);
			const pt::op_t<pt::idx_t> &gt_col_idx =
				g.get_ref_gt_col_idx(sl.ref_idx);
			this->col_to_allele_count_[gt_col_idx.first]++;
			pt::idx_t alt_col_idx =
				this->alt_at_to_unique_alt_idx_.at(i);
			this->genotype_cols_[gt_col_idx.first]
					    [gt_col_idx.second] =
				std::to_string(alt_col_idx);
			// Count this allele occurrence
			this->allele_counts_[alt_col_idx]++;
		}

		// Count reference alleles (genotypes that remained "0")
		for (const auto &col : this->genotype_cols_) {
			for (const auto &gt : col) {
				if (gt == "0") {
					this->allele_counts_[0]++;
				}
			}
		}
	}

public:
	// --------------
	// constructor(s)
	// --------------

	VcfRec(pt::id_t ref_id, pt::idx_t pos, std::string id,
	       pga::allele_slice_t ref_at, pt::idx_t height, var_type_e var_typ,
	       bool is_tangled, pt::idx_t ref_at_ref_count,
	       std::vector<std::vector<std::string>> &&genotype_cols)
	    : ref_id_(ref_id), pos_(pos), id_(id), ats_({ref_at}),
	      height_(height), var_type_(var_typ), is_tangled_(is_tangled),
	      ref_count(std::vector<pt::idx_t>(ref_at_ref_count)),
	      genotype_cols_(std::move(genotype_cols))
	{}

	// because of allele slice let's make these operators/methods explicit
	VcfRec(const VcfRec &) = delete;
	VcfRec &operator=(const VcfRec &) = delete;
	VcfRec(VcfRec &&) noexcept = default;
	VcfRec &
	operator=(VcfRec &&) = delete; // or = default if you later allow it
	~VcfRec() = default;

	// ---------
	// getter(s)
	// ---------
	pt::id_t get_ref_id() const
	{
		return this->ref_id_;
	}

	pt::idx_t get_pos() const
	{
		return this->pos_;
	}

	const std::string &get_id() const
	{
		return this->id_;
	}

	const std::string &get_qual() const
	{
		return this->qual;
	}

	const std::string &get_filter() const
	{
		return this->filter;
	}

	const std::string &get_format() const
	{
		return this->format;
	}

	const pga::allele_slice_t &get_ref_at() const
	{
		return this->ats_.at(0);
	}

	const pga::allele_slice_t &get_slice(pt::idx_t idx) const
	{
		return this->ats_.at(idx);
	}

	// TODO: if C++ 20 use span
	std::vector<pt::idx_t> get_unique_alt_idxs() const
	{
		std::vector<pt::idx_t> alts_;
		alts_.reserve(this->unique_alts_.size() - 1);
		for (pt::idx_t i = 1; i < this->unique_alts_.size(); ++i) {
			alts_.emplace_back(this->unique_alts_[i]);
		}
		return alts_;
	}

	pt::idx_t get_height() const
	{
		return this->height_;
	}

	var_type_e get_var_type() const
	{
		return this->var_type_;
	}

	bool is_tangled() const
	{
		return this->is_tangled_;
	}

	/* info */
	std::string get_ac() const
	{
		std::string ac_str;
		std::string sep = "";
		// AC field shows counts for alternate alleles only (skip ref at
		// index 0)
		for (pt::idx_t i = 1; i < this->allele_counts_.size(); ++i) {
			ac_str += sep;
			ac_str += std::to_string(this->allele_counts_[i]);
			sep = ",";
		}
		return ac_str;
	}

	std::string get_af() const
	{
		std::string af_str;
		std::string sep = "";

		const pt::idx_t AN = this->get_an();

		if (AN == 0) {
			ERR("AN should never be 0");
			std::exit(EXIT_FAILURE);
		}

		// AF field shows frequencies for alternate alleles only (skip
		// ref at index 0)
		for (pt::idx_t i = 1; i < this->allele_counts_.size(); ++i) {
			af_str += sep;
			double v =
				static_cast<double>(this->allele_counts_[i]) /
				static_cast<double>(AN);
			af_str += fmt::format("{:.1f}", v);
			sep = ",";
		}
		return af_str;
	}

	pt::idx_t get_an() const
	{
		// Total number of alleles = sum of all allele counts
		return std::reduce(this->allele_counts_.begin(),
				   this->allele_counts_.end());
	}

	pt::idx_t get_ns() const
	{
		// sum of values in col_to_allele_count_
		pt::idx_t ns{0};
		for (const auto &[_, allele_count] :
		     this->col_to_allele_count_) {
			ns += allele_count;
		}
		return ns;
	}

	std::string get_genotype_fields() const
	{
		// concatenate each column with tab
		std::ostringstream os;
		for (pt::idx_t i{}; i < this->genotype_cols_.size(); ++i) {
			if (i != 0) {
				os << "\t";
			}
			os << pu::concat_with(this->genotype_cols_[i], '|');
		}
		return os.str();
	}

	// ---------
	// setter(s)
	// ---------
	pt::idx_t append_alt_at(const pga::allele_slice_t &alt_at,
				pt::idx_t alt_at_ref_count)
	{
		this->ats_.emplace_back(alt_at);
		this->ref_count.emplace_back(alt_at_ref_count);
		return this->ats_.size() - 1;
	}

	void gen_rec_data_lookups(const bd::VG &g)
	{
		this->comp_unique_alt_ats();
		this->set_genotype_data(g);
	}
};

// TODO: [c] find a better name
class VcfRecIdx
{
	std::map<pt::idx_t, std::vector<VcfRec>> vcf_recs_;

public:
	// --------------
	// constructor(s)
	// --------------

	VcfRecIdx() : vcf_recs_()
	{}

	// ---------
	// getter(s)
	// ---------
	std::map<pt::idx_t, std::vector<VcfRec>> &get_recs_mut()
	{
		return this->vcf_recs_;
	}

	// create *on purpose* if missing
	std::vector<VcfRec> &ensure_recs_mut(pt::idx_t ref_id)
	{
		// inserts empty vector if missing
		auto [it, inserted] = vcf_recs_.try_emplace(ref_id);
		return it->second;
	}

	// ---------
	// setter(s)
	// ---------

	void
	ensure_append_recs(std::map<pt::idx_t, std::vector<VcfRec>> &&new_recs)
	{
		for (auto &&[ref_id, recs] : new_recs) {
			auto &d = this->vcf_recs_[ref_id];
			d.reserve(d.size() + recs.size());
			for (auto &r : recs) {
				d.emplace_back(std::move(r));
			}
		}
	}
};

VcfRecIdx gen_vcf_records(const bd::VG &g,
			  const std::vector<pga::Exp> &ref_walks,
			  const std::set<pt::id_t> &to_call_ref_ids);

} // namespace povu::genomics::vcf

#endif // POVU_GENOMICS_VCF_HPP
