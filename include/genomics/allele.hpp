#ifndef POVU_GENOMICS_ALLELE_HPP
#define POVU_GENOMICS_ALLELE_HPP

#include <algorithm>
#include <deque>
#include <functional> // for std::reference_wrapper
#include <future>
#include <optional>
#include <set>
#include <stack>
#include <string>
#include <sys/types.h>
#include <utility>
#include <vector>

#include "../common/compat.hpp"
#include "../common/types/core.hpp"
#include "../common/types/genomics.hpp"
#include "../common/types/graph.hpp"
#include "../common/types/pvst.hpp"
#include "../common/log.hpp"
#include "../graph/bidirected.hpp"

namespace povu::genomics::allele {
inline constexpr std::string_view MODULE = "povu::genomics::allele";

namespace pvt = povu::types::genomics;

class WalkRefIdx {
  std::set<pt::idx_t> refs_;

  // --------------
  // constructor(s)
  // --------------
  WalkRefIdx() {}

public:
  // --------------
  // factory method(s)
  // --------------
  static WalkRefIdx from_walk(const bd::VG &g, const pvt::walk_t &w) {
    WalkRefIdx wr_idx;
    for (const auto &[v_id, _] : w) {
      const bd::VtxRefIdx &vr_idx = g.get_vtx_ref_idx(v_id);
      wr_idx.refs_.insert(vr_idx.get_ref_ids().begin(), vr_idx.get_ref_ids().end());
    }

    return wr_idx;
  }

  // --------------
  // getter(s)
  // --------------
  const std::set<pt::idx_t> &get_ref_ids() const { return this->refs_; }


  pt::idx_t loop_count(const bd::VtxRefIdx &vri, pt::id_t ref_id) const {
    if (!pv_cmp::contains(this->refs_, ref_id)) {
      return 0;
    }

    if (auto opt_loci = vri.get_ref_loci(ref_id); opt_loci) {
      return opt_loci->get().size(); // the set is expected to be non-empty
    }

    return 0;
  }

  std::optional<pt::idx_t> get_min_locus(const bd::VtxRefIdx &vri, pt::idx_t ref_id) const {
    if (!pv_cmp::contains(this->refs_, ref_id)) {
      return std::nullopt;
    }

    if (auto opt_loci = vri.get_ref_loci(ref_id); opt_loci) {
      return *opt_loci->get().begin(); // the set is expected to be non-empty
    }

    return std::nullopt;
  }
};

void comp_itineraries(const bd::VG &g, const pvt::walk_t &w, pt::idx_t w_idx, pvt::Exp &rw);

} // namespace povu::genomics::alele

#endif // POVU_GENOMICS_ALLELE_HPP
