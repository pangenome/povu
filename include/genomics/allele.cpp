#include "./allele.hpp"
#include <cstdlib>
#include <map>
#include <optional>
#include <sys/types.h>
#include <thread>
#include <unordered_set>
#include <vector>
#include <valarray>


namespace povu::genomics::allele {


bool is_contained(std::vector<pt::slice_t> ref_slices, pt::slice_t ref_slice) {
  for (const pt::slice_t &s : ref_slices) {

    // check if is prefix
    if (s.start == ref_slice.start && s.len > ref_slice.len) {
      return true;
    }

    // check if is suffix
    if (s.start < ref_slice.start && s.start + s.len == ref_slice.start + ref_slice.len) {
      return true;
    }

    // check if is contained
    if (s.start > ref_slice.start && s.start + s.len < ref_slice.start + ref_slice.len) {
      return true;
    }

    // if (s.start <= ref_slice.start && s.start + s.len >= ref_slice.start + ref_slice.len) {
    //   return false;
    // }
  }

  return false;
}

pt::idx_t is_valid(const pgt::ref_walk_t &ref_w, pt::idx_t ref_w_start_idx,
                   const pgt::walk_t &w, pt::idx_t w_start_idx, pt::idx_t len) {
  pt::idx_t valid_len{0};

  for (pt::idx_t i{}; i < len; i++) {
    pt::idx_t ref_w_idx = ref_w_start_idx + i;
    pt::idx_t graph_w_idx = w_start_idx + i;

    auto [ref_v_id, ref_o, _] = ref_w[ref_w_idx];
    auto [w_v_id, w_o] = w[graph_w_idx];

    if (ref_v_id != w_v_id || ref_o != w_o) {
      break;
    }

    valid_len++;
  }

  return valid_len;
}

bool comp_overlays(const bd::VG &g, const pgt::walk_t &w, pt::idx_t w_idx,
                   std::map<pt::id_t, itn_t> &ref_map,
                   std::map<pt::idx_t, std::set<pt::idx_t>> &walk_to_refs) {

  bool is_tangled {false};
  const pt::idx_t WALK_LEN = w.size();

  for (pt::idx_t ref_idx{}; ref_idx < g.get_ref_count(); ref_idx++) {
    const pgt::ref_walk_t &ref_w = g.get_ref_vec(ref_idx);

    itn_t ref_itns;

    std::vector<pt::slice_t> walk_slices;

    for (pt::idx_t walk_step_idx{}; walk_step_idx < WALK_LEN; ++walk_step_idx) {
      pt::idx_t slice_len = WALK_LEN - walk_step_idx;
      const pgt::step_t &step = w[walk_step_idx];
      auto [v_id, o] = step;
      pt::idx_t v_idx = g.v_id_to_idx(v_id);
      const std::vector<pt::idx_t> &vtx_ref_idxs = g.get_vertex_ref_idxs(v_idx, ref_idx);

      std::pair<const pgt::walk_t&, pt::slice_t> walk_slice = {w, {walk_step_idx, slice_len}};
      for (pt::idx_t ref_step_idx : vtx_ref_idxs) {

        // Assumption:
        // in a valid GFA file,
        // if a ref starts within a walk then it has to start at the beginning
        // of the walk
        if (ref_step_idx != 0 && walk_step_idx > 0) {
          continue;
        }

        std::pair<const pgt::ref_walk_t &, pt::slice_t> ref_slice = {ref_w, {ref_step_idx, slice_len}};

        pt::idx_t valid_len = is_valid(ref_w, ref_step_idx, w, walk_step_idx, slice_len);

        if (valid_len < 2) {
          continue;
        }

        if (!is_contained(walk_slices, {walk_step_idx, valid_len})) {
          walk_slices.push_back({walk_step_idx, valid_len});
          ref_itns.append_at({&w, w_idx, walk_step_idx, &ref_w, ref_idx, ref_step_idx, valid_len});
          walk_to_refs[w_idx].insert(ref_idx);

          if (ref_itns.at_count() > 1) {
            is_tangled = true;
          }
        }
      }
    }

    if (ref_itns.at_count() > 0) {
      ref_map.emplace(ref_idx, std::move(ref_itns));
    }

  }
  return is_tangled;
}



void comp_itineraries(const bd::VG &g, Exp &exp) {

  if (exp.get_rov() == nullptr) {
    ERR("RoV pointer is null");
    std::exit(EXIT_FAILURE);
  }

  //bool dbg = exp.id() == ">3645>3647" ? true : false;

  const std::vector<pgt::walk_t> &walks = exp.get_rov()->get_walks();
  std::map<pt::id_t, itn_t> &ref_map = exp.get_ref_itns_mut();
  std::map<pt::idx_t, std::set<pt::idx_t>> &walk_to_refs = exp.get_walk_to_ref_idxs_mut();

  for (pt::idx_t w_idx = 0; w_idx < walks.size(); ++w_idx) {

    // if (dbg) {
    //   INFO("{}", pgt::to_string(walks[w_idx]));
    // }

    bool is_tangled = comp_overlays(g, walks.at(w_idx), w_idx, ref_map, walk_to_refs);
    if (is_tangled) {
      exp.set_tangled(true);
    }

#ifdef DEBUG
    auto refs = exp.get_ref_idxs_for_walk(w_idx); // ensure walk pointers are valid
    for (auto r : refs) {
      const itn_t &r_itn = exp.get_itn(r);
      for (pt::idx_t at_idx = 0; at_idx < r_itn.at_count(); ++at_idx) {
        const allele_slice_t &as = r_itn.get_at(at_idx);
        auto w_idx = as.walk_idx;
        const pgt::walk_t &w = exp.get_rov()->get_walks().at(w_idx);
        if (as.walk != &w) {
          ERR("Inconsistent walk pointer in allele_slice, {}", exp.id());
          std::exit(EXIT_FAILURE);
        }

        pt::idx_t step_idx = as.walk_start_idx;
        pt::idx_t end = as.walk_start_idx + as.len;
        for (step_idx; step_idx < end; ++step_idx) {
          auto _ = as.get_step(step_idx);
        }
      }
    }
#endif

  }

  return;
}

} // namespace povu::genomics::allele
