#ifndef OZ_REFS_HPP
#define OZ_REFS_HPP

#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <string_view>

#include <liteseq/gfa.h>
#include <log/log.h>
#include <quilt/constants.hpp>	 // for
#include <quilt/graph_types.hpp> // for v_end_e, side_n_id_t, side_n_idx_t
#include <quilt/shim.hpp>	 // for contains
#include <quilt/types.hpp>	 // for qt
#include <quilt/utils.hpp>	 // for pu, TwoWayMap

namespace oza::refs
{
namespace lq = liteseq;

enum class ref_format_e : qt::u8 {
	PANSN = 1,
	UNDEFINED = 0,
};
std::string to_string(ref_format_e r);

char lq_strand_to_char(liteseq::strand s);

ptg::or_e lq_strand_to_pv_or(liteseq::strand s);

/**
 * zero indexed
 * hap col 0 is the first haplotype column
 * phase col 0 is the first phase column
 */
struct gt_col_meta {
	qt::u32 hap_col;   // haplotype col idx
	qt::u32 phase_col; // phase col idx or ploidy
};

struct ploidy_meta {
private:
	// hap ids of the sample
	// idx 0 is the first hap_id, idx 1 is the second hap id, etc.
	std::vector<qt::u32> hap_ids;

public:
	// -----------
	// constructor
	// -----------
	ploidy_meta() = default;

	// -------
	// getters
	// -------
	[[nodiscard]]
	qt::u32 ploidy() const
	{
		return static_cast<qt::u32>(this->hap_ids.size());
	}

	qt::u32 get_ploidy_idx(qt::u32 ploidy_id) const
	{
		for (qt::u32 i = 0; i < this->hap_ids.size(); i++)
			if (this->hap_ids[i] == ploidy_id)
				return i;

		std::string contents = pu::concat_with(this->hap_ids, ',');
		std::string err = qs::format(
			"Hap id {} not found in ploidy metadata. Contains: {}",
			ploidy_id, contents);
		log_fatal("%s", err.c_str());
		// PL_ERR("Hap id {} not found in ploidy metadata. Contains:
		// {}",
		//        ploidy_id, contents);
		std::exit(EXIT_FAILURE);
	}

	[[nodiscard]]
	qt::u32 get_hap_id(qt::u32 ploidy_idx) const
	{
		if (ploidy_idx >= this->hap_ids.size()) {
			std::string err = qs::format(
				"Hap idx {} out of bounds for ploidy {}",
				ploidy_idx, this->hap_ids.size());
			log_fatal("%s", err.c_str());
			// PL_ERR("Hap idx {} out of bounds for ploidy {}",
			//        ploidy_idx, this->hap_ids.size());
			std::exit(EXIT_FAILURE);
		}

		return this->hap_ids.at(ploidy_idx);
	}

	// -------
	// setters
	// -------
	void add_hap_id(qt::u32 hap_id)
	{
		// adds before the first greater than hap id to keep sorted
		auto it = std::lower_bound(this->hap_ids.begin(),
					   this->hap_ids.end(), hap_id);

		if (it != this->hap_ids.end() && *it == hap_id)
			return; // Duplicate found, exit early

		this->hap_ids.insert(it, hap_id);
	}
};

class Ref
{
	const lq::ref *ref_ptr_; // pointer to the original liteseq ref

	// make default constructor private
	Ref() = default;

public:
	// --------------
	// constructor(s)
	// --------------
	static Ref from_lq_ref(const lq::ref *r)
	{
		Ref ref;
		ref.ref_ptr_ = r;
		return ref;
	}

	// ---------
	// getter(s)
	// ---------
	[[nodiscard]]
	std::string tag() const
	{
		const char *t = lq::get_tag(this->ref_ptr_);
		if (t == nullptr)
			std::cerr << "Failed to get tag for ref\n";
		return t;
	}

	[[nodiscard]]
	qt::id_t get_hap_id() const
	{
		return lq::get_hap_id(this->ref_ptr_);
	}

	[[nodiscard]]
	ref_format_e get_format() const
	{
		lq::ref_id_type rt = lq::get_ref_id_type(this->ref_ptr_);
		if (rt == lq::ref_id_type::REF_ID_PANSN) {
			return ref_format_e::PANSN;
		}
		else {
			return ref_format_e::UNDEFINED;
		}
	}

	[[nodiscard]]
	std::string get_sample_name() const
	{
		return lq::get_sample_name(this->ref_ptr_);
	}

	[[nodiscard]]
	qt::id_t get_length() const
	{
		return lq::get_hap_len(this->ref_ptr_);
	}
};

class Refs
{
	lq::ref **refs_ptr_ptr;

	qt::idx_t ref_count_;

	// TODO we are storing this pointer twice. Let's not do that.
	std::vector<Ref> refs_;

	// sn2pm: sample name to ploidy metadata
	std::map<std::string, std::optional<ploidy_meta>> sn2pm;

	// sample name to hap idxs
	std::map<std::string, std::set<qt::u32>> sn2h_idxs;

	/*
	  -------------
	  genotype data
	  -------------
	*/

	// TODO: in C++20 or greater, use a constexpr instead of const
	inline static const std::string BLANK_GT_VALUE = ".";

	std::vector<std::vector<std::string>> blank_genotype_cols;

	// map from hap id to genotype column indices
	std::map<qt::idx_t, gt_col_meta> ref_id_to_col_idx;

	// TODO: [remove] we don't need to store this
	// the size is the number of columns in the genotype data
	// string contains the sample name or label
	std::vector<std::string> genotype_col_names;

public:
	// ---------
	// getter(s)
	// ---------
	[[nodiscard]]
	qt::id_t ref_count() const
	{
		return ref_count_;
	}

	[[nodiscard]]
	const lq::ref *get_lq_ref_ptr(qt::idx_t ref_idx) const
	{
		return refs_ptr_ptr[ref_idx];
	}

	[[nodiscard]]
	const Ref &get_lq_ref(qt::id_t ref_id) const
	{
		return this->refs_.at(ref_id);
	}

	[[nodiscard]]
	std::string get_sample_name(qt::id_t ref_id) const
	{
		const Ref &r = this->refs_.at(ref_id);
		std::string sn = r.get_sample_name();
		return sn;
	}

	[[nodiscard]]
	qt::u32 get_ploidy(const std::string &sample_name) const
	{
		auto it = this->sn2pm.find(sample_name);
		if (it == this->sn2pm.end()) {
			std::string err = qs::format("Sample name {} not found",
						     sample_name);
			log_fatal("%s", err.c_str());
			// PL_ERR("Sample name {} not found", sample_name);
			std::exit(EXIT_FAILURE);
		}

		// TODO: handle undefined ploidy more gracefully
		// this is for non PANSN refs
		const std::optional<ploidy_meta> &pm = it->second;
		if (!pm.has_value()) {
			return pc::INVALID_IDX;
			// PL_ERR("Sample name {} has undefined ploidy",
			//        sample_name);
			// std::exit(EXIT_FAILURE);
		}

		return pm.value().ploidy();
	}

	[[nodiscard]]
	qt::u32 get_ploidy_id(const std::string &sample_name,
			      qt::u32 ploidy_idx) const
	{
		auto it = this->sn2pm.find(sample_name);
		if (it == this->sn2pm.end()) {
			std::string err = qs::format("Sample name {} not found",
						     sample_name);
			log_fatal("%s", err.c_str());
			std::exit(EXIT_FAILURE);
		}

		const std::optional<ploidy_meta> &pm = it->second;
		if (!pm.has_value()) {
			return pc::INVALID_ID;
			// PL_ERR("Sample name {} has undefined ploidy at {}",
			//        sample_name, ploidy_idx);
			// std::exit(EXIT_FAILURE);
		}

		return pm.value().get_hap_id(ploidy_idx);
	}

	[[nodiscard]]
	std::string get_tag(qt::id_t ref_id) const
	{
		const Ref &r = this->refs_.at(ref_id);
		return r.tag();
	}

	[[nodiscard]]
	std::optional<qt::id_t> get_ref_id(std::string_view tag) const
	{
		for (qt::id_t ref_id{}; ref_id < this->ref_count(); ref_id++) {
			const char *t = lq::get_tag(this->refs_ptr_ptr[ref_id]);
			if (tag == t)
				return ref_id;
		}

		return std::nullopt;
	}

	// the sample name could also be referred to as a prefix
	[[nodiscard]]
	std::set<qt::id_t>
	get_refs_in_sample(std::string_view sample_name) const
	{
		std::string sn{sample_name};
		std::set<qt::u32> combined;
		auto it = this->sn2h_idxs.lower_bound(sn);
		while (it != this->sn2h_idxs.end() &&
		       pu::is_prefix(sn, it->first)) {
			combined.insert(it->second.begin(), it->second.end());
			++it;
		}

		return combined;
	}

	[[nodiscard]]
	std::set<qt::id_t> get_shared_samples(qt::id_t ref_id) const
	{
		const lq::ref *r = this->refs_ptr_ptr[ref_id];
		const char *sample_name = lq::get_sample_name(r);
		return this->get_refs_in_sample(std::string_view(sample_name));
	}

	[[nodiscard]]
	const std::vector<std::string> &get_genotype_col_names() const
	{
		return this->genotype_col_names;
	}

	[[nodiscard]]
	const std::vector<std::vector<std::string>> &
	get_blank_genotype_cols() const
	{
		return this->blank_genotype_cols;
	}

	[[nodiscard]]
	const gt_col_meta &get_gt_col_idx(qt::id_t ref_id) const
	{
		return this->ref_id_to_col_idx.at(ref_id);
	}

	// ---------
	// setter(s)
	// ---------
	/**
	 * @brief adds ref names and ref metadata
	 *
	 */
	void set_refs_meta(lq::ref **refs, qt::u32 ref_count)
	{
		this->ref_count_ = ref_count;
		this->refs_ptr_ptr = refs;
		for (qt::idx_t h_idx{}; h_idx < ref_count; h_idx++) {
			Ref ref = Ref::from_lq_ref(refs_ptr_ptr[h_idx]);
			this->refs_.emplace_back(ref);

			std::string sn = ref.get_sample_name();

			if (!qs::contains(this->sn2h_idxs, sn))
				this->genotype_col_names.emplace_back(sn);

			std::set<qt::u32> &sample_hap_idxs =
				this->sn2h_idxs[sn];
			sample_hap_idxs.insert(h_idx);
		}

		// set ploidy metadata. ploidy meta only makes sense for PanSN
		for (const auto &[sn, h_idxs] : this->sn2h_idxs) {
			std::optional<ploidy_meta> &opt_pm =
				this->sn2pm[sn]; // Create nullopt if none

			for (qt::u32 h_idx : h_idxs) {
				const Ref &r = this->get_lq_ref(h_idx);

				if (r.get_format() != ref_format_e::PANSN)
					continue;

				qt::u32 ploidy_id = r.get_hap_id();
				if (!opt_pm)
					opt_pm = ploidy_meta{};

				opt_pm->add_hap_id(ploidy_id);
			}
		}

		// compute and save genotype fields for each hap
		qt::u32 N = this->genotype_col_names.size();
		for (qt::u32 hap_col{}; hap_col < N; hap_col++) {
			const std::string &sn = genotype_col_names[hap_col];
			const std::set<qt::u32> &sample_haps = sn2h_idxs[sn];
			const std::optional<ploidy_meta> &opt_pm = sn2pm[sn];

			// ---
			// pre compute the blank gt cols
			// if it has ploidy meta set ploidy to what is computed
			// else we assume that the number of sample haps is the
			// ploidy
			// ----

			qt::u32 ploidy = (opt_pm) ? this->get_ploidy(sn)
						  : sample_haps.size();
			this->blank_genotype_cols.emplace_back(ploidy,
							       BLANK_GT_VALUE);

			// ---
			//
			// ----
			qt::u32 phase_col_ctr =
				0; // used when there is no phase info
			for (qt::u32 h_idx : sample_haps) {
				const Ref &r = this->get_lq_ref(h_idx);
				ref_format_e ref_fmt = r.get_format();

				if (ref_fmt == ref_format_e::UNDEFINED) {
					this->ref_id_to_col_idx[h_idx] = {
						hap_col, phase_col_ctr++};
				}
				else if (ref_fmt == ref_format_e::PANSN) {
					qt::u32 phase_id = r.get_hap_id();
					qt::u32 phase_col =
						opt_pm->get_ploidy_idx(
							phase_id);

					this->ref_id_to_col_idx[h_idx] = {
						hap_col, phase_col};
				}
			}
		}
	}
};

}; // namespace oza::refs

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pr = oza::refs;

#endif // OZ_REFS_HPP
