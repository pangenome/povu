#ifndef POVU_GRAPH_REF_HPP
#define POVU_GRAPH_REF_HPP

#include <map>
#include <optional>
#include <set>
#include <string_view>


#include "../common/core.hpp"
#include "./types.hpp"


namespace povu::graph::ref {
inline constexpr std::string_view MODULE = "povu::graph::ref";

namespace ptg = povu::types::graph;

struct ref_step_t {
  ptg::or_e strand_;
  pt::idx_t locus_;

  ref_step_t(ptg::or_e strand, pt::idx_t locus) : strand_(strand), locus_(locus) {}
};
bool operator==(const ref_step_t &lhs, const ref_step_t &rhs);
bool operator<(const ref_step_t &lhs, const ref_step_t &rhs);


class VtxRefData {
  std::set<pt::id_t> ref_ids_;
  // key is a  ref id the value is the set of steps the ref makes in the vertex
  // sorted by locus
  std::map<pt::id_t, std::set<ref_step_t>> ref_id_to_data_;

public:
  // --------------
  // constructor(s)
  // --------------
  VtxRefData() = default;

  // ---------
  // getter(s)
  // ---------

  bool has_ref(pt::id_t ref_id) const {
    return pv_cmp::contains(this->ref_id_to_data_, ref_id);
  }

  bool has_address(pt::id_t ref_id, ref_step_t a) const {
    if (!this->has_ref(ref_id)) {
      return false;
    }
    const std::set<ref_step_t> &data = this->ref_id_to_data_.at(ref_id);
    return pv_cmp::contains(data, a);
  }

  const std::set<ref_step_t> &get_ref_data(pt::id_t ref_id) const {
    return this->ref_id_to_data_.at(ref_id);
  }

  const std::set<pt::id_t> &get_ref_ids() const {
    return this->ref_ids_;
  }

  pt::idx_t loop_count(pt::id_t ref_id) const {
    if (!this->has_ref(ref_id)) {
      return 0;
    }
    return this->get_ref_data(ref_id).size();
  }

  // TODO: should we return the entire step?
  // works with 0 because we assume indexes are 1 indexed
  // the min must be > threshold
  std::optional<pt::idx_t>
  get_min_locus(pt::id_t ref_id, std::optional<pt::idx_t> opt_start_after) const {
    if (!pv_cmp::contains(this->ref_ids_, ref_id)) {
      return std::nullopt;
    }

    pt::idx_t threshold = opt_start_after.value_or(0);
    for (const ref_step_t &a : this->ref_id_to_data_.at(ref_id)) {
      if (a.locus_ > threshold) {
        return a.locus_;
      }
    }
    return std::nullopt;
  }

  // ---------
  // setter(s)
  // ---------

  void add_ref_data(pt::id_t ref_id, ptg::or_e strand, pt::idx_t locus) {
    this->ref_id_to_data_[ref_id].insert({strand, locus});
    this->ref_ids_.insert(ref_id);
  }
};

} // namespace povu::graph::ref

// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pgr = povu::graph::ref;

#endif // POVU_GRAPH_REF_HPP
