#ifndef POVU_REFS_HPP
#define POVU_REFS_HPP

#include <liteseq/gfa.h>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <string_view>

#include "povu/common/core.hpp"
#include "povu/common/log.hpp"
#include "povu/common/utils.hpp"

namespace povu::refs
{
inline constexpr std::string_view MODULE = "povu::refs";
namespace lq = liteseq;

enum class ref_format_e {
	PANSN = 1,
	UNDEFINED = 0,
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
	std::string tag() const
	{
		return lq::get_tag(this->ref_ptr_);
	}

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

	const std::string &get_sample_name() const
	{
		return lq::get_sample_name(this->ref_ptr_);
	}

	pt::id_t get_length() const
	{
		return lq::get_hap_len(this->ref_ptr_);
	}
};

class Refs
{
	lq::ref **refs;

	pt::idx_t ref_count_;

	// TODO we are storing this pointer twice. Let's not do that.
	std::vector<Ref> refs_;

	/* genotype data */

	// TODO: in C++20 or greater, use a constexpr instead of const
	inline static const std::string BLANK_GT_VALUE = ".";

	std::vector<std::vector<std::string>> blank_genotype_cols;

	// std::map<pt::id_t, pt::idx_t> ref_id_to_col_idx;
	//  ref id to two dimensional col idx
	std::map<pt::idx_t, pt::op_t<pt::idx_t>> ref_id_to_col_idx;

	// TODO: [remove] we don't need to store this
	// the size is the number of columns in the genotype data
	// string contains the sample name or label
	std::vector<std::string> genotype_col_names;

public:
	// ---------
	// getter(s)
	// ---------

	pt::id_t ref_count() const
	{
		return ref_count_;
	}

	const lq::ref *get_lq_ref_ptr(pt::idx_t ref_idx) const
	{
		return refs[ref_idx];
	}

	const Ref &get_lq_ref(pt::id_t ref_id) const
	{
		return this->refs_.at(ref_id);
	}

	const std::string &get_sample_name(pt::id_t ref_id) const
	{
		return lq::get_sample_name(this->refs[ref_id]);
	}

	std::optional<pt::id_t> get_ref_id(std::string_view tag) const
	{
		for (pt::id_t ref_id{}; ref_id < this->ref_count(); ref_id++) {
			const char *t = lq::get_tag(this->refs[ref_id]);
			if (tag == t)
				return ref_id;
		}

		return std::nullopt;
	}

	// the sample name could also be referred to as a prefix
	std::set<pt::id_t>
	get_refs_in_sample(std::string_view sample_name) const
	{
		std::set<pt::id_t> in_sample;
		for (pt::id_t ref_id{}; ref_id < this->ref_count(); ref_id++) {
			const lq::ref *r = this->refs[ref_id];
			const char *r_sn = lq::get_sample_name(r);
			lq::ref_id_type rt = lq::get_ref_id_type(r);
			if (rt == lq::ref_id_type::REF_ID_PANSN) {
				if (r_sn == sample_name)
					in_sample.insert(ref_id);
			}
			else if (rt == lq::ref_id_type::REF_ID_RAW) {
				if (pu::is_prefix(sample_name, r_sn))
					in_sample.insert(ref_id);
			}
		}

		return in_sample;
	}

	std::set<pt::id_t> get_shared_samples(pt::id_t ref_id) const
	{
		const lq::ref *r = this->refs[ref_id];
		const char *sample_name = lq::get_sample_name(r);
		return this->get_refs_in_sample(sample_name);
	}

	const std::vector<std::string> &get_genotype_col_names() const
	{
		return this->genotype_col_names;
	}

	const std::vector<std::vector<std::string>> &
	get_blank_genotype_cols() const
	{
		return this->blank_genotype_cols;
	}

	const pt::op_t<pt::idx_t> &get_ref_gt_col_idx(pt::id_t ref_id) const
	{
		return this->ref_id_to_col_idx.at(ref_id);
	}

	// ---------
	// setter(s)
	// ---------
	void add_all_refs(lq::ref **refs, pt::idx_t ref_count)
	{
		this->ref_count_ = ref_count;
		this->refs = refs;
		for (pt::idx_t ref_idx{}; ref_idx < ref_count; ++ref_idx) {
			Ref ref = Ref::from_lq_ref(refs[ref_idx]);
			this->refs_.emplace_back(std::move(ref));
		}
	}

	void gen_genotype_metadata()
	{
		std::set<pt::id_t> handled;

		for (pt::id_t ref_id = 0; ref_id < this->ref_count();
		     ++ref_id) {
			if (pv_cmp::contains(handled, ref_id)) {
				continue;
			}

			handled.insert(ref_id);

			const std::set<pt::id_t> &sample_refs =
				this->get_shared_samples(ref_id);

			const lq::ref *r = this->refs[ref_id];
			const char *col_name = lq::get_sample_name(r);
			// std::string col_name =
			// this->get_ref(ref_id).get_sample_name();

			if (sample_refs.empty()) {

				ERR("No sample names found for ref_id {}",
				    ref_id);
				std::exit(EXIT_FAILURE);
			}
			else if (sample_refs.size() == 1) {
				this->blank_genotype_cols.push_back(
					std::vector<std::string>{
						BLANK_GT_VALUE});
				pt::op_t<pt::idx_t> x{
					static_cast<pt::idx_t>(
						this->genotype_col_names
							.size()),
					0};
				this->ref_id_to_col_idx[ref_id] = x;
				this->genotype_col_names.push_back(col_name);
			}
			else if (sample_refs.size() > 1) {
				pt::idx_t col_idx =
					this->genotype_col_names.size();
				pt::idx_t col_col_idx{};

				this->blank_genotype_cols.push_back(
					std::vector<std::string>(
						sample_refs.size(),
						BLANK_GT_VALUE));

				for (pt::id_t ref_id_ : sample_refs) {
					pt::op_t<pt::idx_t> x{col_idx,
							      col_col_idx++};
					this->ref_id_to_col_idx[ref_id_] = x;
					handled.insert(ref_id_);
				}
				this->genotype_col_names.push_back(col_name);
			}
		}
	}
};

}; // namespace povu::refs

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pr = povu::refs;

#endif // POVU_REFS_HPP
