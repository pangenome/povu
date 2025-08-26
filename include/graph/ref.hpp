#ifndef POVU_GRAPH_REF_HPP
#define POVU_GRAPH_REF_HPP

#include <map>
#include <optional>
#include <set>
#include <string>
#include <vector>


#include "../common/core.hpp"
#include "../common/utils.hpp"


namespace povu::graph::ref {


struct pan_sn {
  std::string sample_name_; // not unique
  pt::id_t haplotype_id_;   // is unique
  std::string contig_name_; // is unique contig or scaffold name
};

std::optional<pan_sn> label_to_pan_sn(const std::string &label, char delim);

// applies to a single line in the P line of a GFA file

class Ref {
  pt::id_t ref_id_;
  std::string label_; // the entire label in the P line

  pt::idx_t len_; // length of the reference/contig/scaffold

  // PanSN spec
  bool has_pansn_data_ = false; // true if the label is in the PanSN format
  pan_sn pansn_data_; // holds the PanSN data if has_pansn_data_ is true
  // std::string sample_name_; // not unique
  // pt::id_t haplotype_id_; // is unique
  // std::string contig_name_; // is unique contig or scaffold name

public:

  // --------------
  // constructor(s)
  // --------------

  static Ref parse_pansn_label(pt::id_t ref_id, const std::string &label, char delim) {
    Ref info;
    info.ref_id_ = ref_id;
    info.label_ = label;

    std::optional<pan_sn> maybe_pansn = label_to_pan_sn(label, delim);
    if (maybe_pansn.has_value()) {
      info.has_pansn_data_ = true;
      info.pansn_data_ = maybe_pansn.value();
    }
    else {
      info.has_pansn_data_ = false;
    }

    return info;
  }

  // ---------
  // getter(s)
  // ---------

  pt::id_t id() const {
    return this->ref_id_;
  }

  bool has_pansn_data() const {
    return this->has_pansn_data_;
  }

  const std::string &get_label() const {
    return this->label_;
  }

  const std::string &get_sample_name() const {
    return this->pansn_data_.sample_name_;
  }

  std::string get_col_name() const {
    if (this->has_pansn_data_) {
      return this->pansn_data_.sample_name_;
    }
    else {
      return this->label_;
    }
  }

  pt::id_t get_haplotype_id() const { return this->pansn_data_.haplotype_id_; }

  const std::string &get_contig_name() const {
    return this->pansn_data_.contig_name_;
  }

  pt::idx_t get_length() const {
    return this->len_;
  }

  // ---------
  // setter(s)
  // ---------

  void set_length(pt::idx_t len) {
    this->len_ = len;
  }
};

class Refs {
  std::vector<Ref> refs_;

  // the ref id is the index in the vector and the value is the set of ref ids
  // which share the sample name with the ref id at that index
  std::vector<std::set<pt::id_t>> ref_ids_;

  // the string is the sample name and the values are the ref ids which share
  // a sample name
  std::map<std::string, std::set<pt::id_t>> sample_to_ref_ids_;

  pu::TwoWayMap<std::string, pt::id_t> label_and_ref_id_;

public:

  // --------------
  // constructor(s)
  // --------------
  Refs() = default;

  // ---------
  // getter(s)
  // ---------

  const Ref &get_ref(pt::id_t ref_id) const {
    return this->refs_.at(ref_id);
  }

  Ref &get_ref_mut(pt::id_t ref_id) {
    return this->refs_.at(ref_id);
  }

  const std::string &get_ref_label(pt::id_t ref_id) const {
    return this->refs_.at(ref_id).get_label();
  }

  pt::id_t get_ref_id(const std::string &ref_label) const {
    return this->label_and_ref_id_.get_value(ref_label);
  }

  pt::id_t ref_id_count() const {
    return this->refs_.size();
  }

  /**
   * Get the set of ref ids that share the same sample name as the ref id
   */
  const std::set<pt::id_t> &get_shared_samples(pt::id_t ref_id) const {
    return this->ref_ids_.at(ref_id);
  }

  // ---------
  // setter(s)
  // ---------

  pt::id_t add_ref(const std::string &label, char delim) {
    pt::id_t ref_id = this->refs_.size();
    Ref ref = Ref::parse_pansn_label(ref_id, label, delim);
    this->refs_.push_back(ref);
    this->label_and_ref_id_.insert(label, ref_id);

    const std::string &sample_name = ref.get_col_name();

    this->sample_to_ref_ids_.try_emplace(sample_name, std::set<pt::id_t>{});
    this->sample_to_ref_ids_[sample_name].insert(ref_id);

    // TODO [B] use already pre computed ref count
    this->ref_ids_.resize(this->refs_.size());

    // update for all refs
    for (pt::id_t r_id : this->sample_to_ref_ids_[sample_name]) {
      this->ref_ids_[r_id] = this->sample_to_ref_ids_[sample_name];
    }

    return ref_id;
  }
};

} // namespace povu::graph::ref

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pgr = povu::graph::ref;

#endif // POVU_GRAPH_REF_HPP
