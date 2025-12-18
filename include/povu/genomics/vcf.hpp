#ifndef POVU_GENOMICS_VCF_HPP
#define POVU_GENOMICS_VCF_HPP

#include <cstdlib> // for exit, EXIT_FAILURE
#include <map>	   // for map, _Rb_tree_iterator, operator!=
// #include <numeric>     // for reduce
#include <set>	       // for set
#include <sstream>     // for basic_ostringstream, basic_ostream
#include <string>      // for basic_string, string, allocator
#include <string_view> // for operator<<, string_view
#include <tuple>       // for tuple, make_tuple
#include <utility>     // for get, move, pair
#include <vector>      // for vector

#include "fmt/core.h" // for format
#include "povu/common/constants.hpp"
#include "povu/common/core.hpp"		  // for pt, idx_t, id_t, op_t
#include "povu/common/log.hpp"		  // for ERR
#include "povu/common/utils.hpp"	  // for concat_with, pu
#include "povu/genomics/allele.hpp"	  // for allele_slice_t, Exp
#include "povu/graph/bidirected.hpp"	  // for VG
#include "povu/overlay/interval_tree.hpp" // for it
#include "povu/variation/rov.hpp"	  // for var_type_e

namespace povu::genomics::vcf
{
inline constexpr std::string_view MODULE = "povu::genomics::vcf";
namespace bd = povu::bidirected;
namespace pgt = povu::types::graph;
namespace pga = povu::genomics::allele;

class VcfRec
{
	pt::id_t ref_id_;	 // chrom TODO: does this still apply with tags?
	pt::idx_t pos_;		 // 1-based step idx
	std::string id_;	 // start and end of a RoV e.g >1>4
	std::string enc_flubble; // enclosing flubble string

	pga::hap_slice ref_slice_;	// the ref allele slice
	std::set<pt::u32> ref_at_haps_; // haps that contain the ref allele

	// the set of alt allele slices
	std::vector<pga::hap_slice> ats_;

	// alt alleles grouped by alt allele idx
	std::vector<std::vector<pga::hap_slice>> alt_slices_;

	pt::u32 ns_{pc::INVALID_IDX};

	inline static const std::string qual = "60";	 // fixed at 60
	inline static const std::string filter = "PASS"; // fixed as pass
	inline static const std::string format = "GT";	 // fixed as pass

	pt::idx_t height_; // height of the pvst node in the tree

	// type of the variant, e.g. del, ins, sub, und
	pvr::var_type_e var_type_;

	// is true when tangling exists, i.e. when a
	// walk traverses an RoV more than once
	bool is_tangled_{false};

	/* info */
	// number of refs in a given walk
	std::vector<pt::idx_t> ref_count;
	std::vector<pt::idx_t> allele_counts_; // count for each allele (ref at
					       // idx 0, alts at idx 1+)
	std::vector<double> af_; // allele frequency for each alt allele
	// pt::idx_t an_ = 0;	 // total number of alleles in called genotypes
	//  pt::idx_t ns_ = 0;       // number of samples with data

	// genotype column idx to allele count
	std::map<pt::idx_t, pt::idx_t> col_to_allele_count_;

	// genotype
	std::vector<std::vector<std::string>> genotype_cols_;

	// lookups
	// each value represents a unique alt allele idx
	std::vector<pt::idx_t> unique_alts_;
	// maps alt allele idx (0-based) to unique alt allele idx (1-based)
	std::vector<pt::idx_t> alt_at_to_unique_alt_idx_;

public:
	// --------------
	// constructor(s)
	// --------------
	VcfRec() = delete;

	VcfRec(pt::u32 ref_h_idx, pt::u32 pos, std::string id,
	       std::string en_flub, pga::hap_slice ref_sl, pt::u32 height,
	       std::set<pt::u32> &&ref_at_haps, pvr::var_type_e var_typ,
	       bool is_tangled)
	    : ref_id_(ref_h_idx), pos_(pos), id_(std::move(id)),
	      enc_flubble(std::move(en_flub)), ref_slice_(ref_sl),
	      ref_at_haps_(ref_at_haps), height_(height), var_type_(var_typ),
	      is_tangled_(is_tangled)
	{}

	// because of allele slice let's make these operators/methods explicit
	VcfRec(const VcfRec &) = delete;
	VcfRec &operator=(const VcfRec &) = delete;
	VcfRec(VcfRec &&) noexcept = default;
	VcfRec &
	operator=(VcfRec &&) = delete; // or = default if you later allow it
	~VcfRec() = default;

	void add_gt_cols(std::vector<std::vector<std::string>> &&genotype_cols)
	{
		this->genotype_cols_ = std::move(genotype_cols);
	}

	// ---------
	// getter(s)
	// ---------
	[[nodiscard]]
	pt::id_t get_ref_id() const

	{
		return this->ref_id_;
	}

	[[nodiscard]]
	pt::idx_t get_pos() const
	{
		return this->pos_;
	}

	[[nodiscard]]
	const std::string &get_id() const
	{
		return this->id_;
	}

	[[nodiscard]]
	const std::string &get_qual() const
	{
		return this->qual;
	}

	[[nodiscard]]
	const std::string &get_filter() const
	{
		return this->filter;
	}

	[[nodiscard]]
	const std::string &get_format() const
	{
		return this->format;
	}

	[[nodiscard]]
	const pga::hap_slice &get_ref_slice() const
	{
		return this->ref_slice_;
	}

	[[nodiscard]]
	std::string get_ref_as_dna_str(const bd::VG &g) const
	{
		return this->get_ref_slice().as_dna_str(g, this->var_type_);
	}

	[[nodiscard]]
	std::string get_alts_as_str() const
	{
		std::string alt_str = "";
		for (pt::u32 i{}; i < this->alt_slices_.size(); ++i) {
			pga::hap_slice alt = this->alt_slices_.at(i).front();
			alt_str += alt.as_str(this->var_type_);
			if (i != this->alt_slices_.size() - 1)
				alt_str += ",";
		}

		return alt_str;
	}

	[[nodiscard]]
	std::string get_alts_as_dna_str(const bd::VG &g) const
	{
		std::string dna_str = "";
		for (pt::u32 i{}; i < this->alt_slices_.size(); ++i) {
			pga::hap_slice alt = this->alt_slices_.at(i).front();
			dna_str += alt.as_dna_str(g, this->var_type_);
			if (i != this->alt_slices_.size() - 1)
				dna_str += ",";
		}

		return dna_str;
	}

	[[nodiscard]]
	std::string get_at() const
	{
		std::string s;
		s += this->get_ref_slice().as_str(this->var_type_);
		s += ",";
		s += this->get_alts_as_str();
		return s;
	}

	[[nodiscard]]
	const pga::hap_slice &get_slice(pt::idx_t idx) const
	{
		return this->ats_.at(idx);
	}

	// TODO: if C++ 20 use span
	[[nodiscard]]
	std::vector<pt::idx_t> get_unique_alt_idxs() const
	{
		std::vector<pt::idx_t> alts_;
		alts_.reserve(this->unique_alts_.size() - 1);
		for (pt::idx_t i = 1; i < this->unique_alts_.size(); ++i)
			alts_.emplace_back(this->unique_alts_[i]);

		return alts_;
	}

	[[nodiscard]]
	pt::idx_t get_height() const
	{
		return this->height_;
	}

	[[nodiscard]]
	std::string get_enc_flubble() const
	{
		return this->enc_flubble;
	}

	[[nodiscard]]
	pvr::var_type_e get_var_type() const
	{
		return this->var_type_;
	}

	[[nodiscard]]
	bool is_tangled() const
	{
		return this->is_tangled_;
	}

	/* info */
	[[nodiscard]]
	std::string get_ac() const
	{
		std::string ac_str;
		std::string sep = "";

		// AC field shows counts for alternate alleles only
		for (pt::u32 i{}; i < this->alt_slices_.size(); i++) {
			ac_str += sep;
			ac_str += std::to_string(this->alt_slices_[i].size());
			sep = ",";
		}
		return ac_str;
	}

	[[nodiscard]]
	std::string get_af() const
	{
		std::string af_str;
		std::string sep = "";

		const pt::idx_t AN = this->get_an();

		if (AN == 0) {
			ERR("AN should never be 0");
			std::exit(EXIT_FAILURE);
		}

		for (pt::u32 i{}; i < this->alt_slices_.size(); i++) {
			af_str += sep;
			double v = static_cast<double>(
					   this->alt_slices_[i].size()) /
				   static_cast<double>(AN);
			af_str += fmt::format("{:.1f}", v);
			sep = ",";
		}

		return af_str;
	}

	[[nodiscard]]
	pt::idx_t get_an() const
	{
		pt::u32 an{};
		// sum the refs
		an += this->ref_at_haps_.size(); // ref allele count

		// sum the alts
		for (pt::u32 i{}; i < this->alt_slices_.size(); i++)
			an += this->alt_slices_[i].size();

		// Total number of alleles = sum of all allele counts
		return an;
	}

	/* Number of samples with data */
	[[nodiscard]]
	pt::idx_t get_ns() const
	{
		return this->ns_;
	}

	void set_ns(pt::u32 ns)
	{
		this->ns_ = ns;
	}

	[[nodiscard]]
	std::string get_genotype_fields() const
	{
		// concatenate each column with tab
		std::ostringstream os;
		for (pt::idx_t i{}; i < this->genotype_cols_.size(); ++i) {
			if (i > 0)
				os << "\t";

			os << pu::concat_with(this->genotype_cols_[i], '|');
		}
		return os.str();
	}

	// ---------
	// setter(s)
	// ---------
	pt::idx_t append_alt_at(const pga::hap_slice &alt_at,
				pt::idx_t alt_at_ref_count)
	{
		this->ats_.emplace_back(alt_at);
		this->ref_count.emplace_back(alt_at_ref_count);
		return this->ats_.size() - 1;
	}

	void add_alt_set(std::vector<pga::hap_slice> alt_set)
	{
		this->alt_slices_.emplace_back(std::move(alt_set));
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

VcfRecIdx gen_vcf_records(const bd::VG &g, const std::vector<pga::trek> &treks,
			  const std::vector<poi::it> &its);

} // namespace povu::genomics::vcf

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pgv = povu::genomics::vcf;

#endif // POVU_GENOMICS_VCF_HPP
