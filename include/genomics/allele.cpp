#include "./allele.hpp"

namespace povu::genomics::allele {


pt::idx_t get_vtx_len(const bd::VG &g, const pvt::step_t &s) {
  const bd::Vertex &v = g.get_vertex_by_id(s.v_id);
  return v.get_label().length();
}

/**
 * Find the minimum locus for a given reference in a walk.
 *
 * The function searches through the walk steps for the smallest locus
 * associated with `ref_id` according to the data in `wri`.
 *
 * Behavior:
 *  - If `opt_start_after` is not provided, the function returns the global
 *    minimum locus across all steps.
 *  - If `opt_start_after` has a value, the function returns the smallest
 *    locus that is strictly greater than `*opt_start_after`. If no such locus
 *    exists, `pc::MAX_IDX` is returned.
 *
 * @param ref_id        The reference ID to search for.
 * @param wri           The WalkRefIdx providing locus information for vertices.
 * @param w             The walk (sequence of steps) to search through.
 * @param opt_start_after Optional threshold value: only loci greater than this
 *                        will be considered. If std::nullopt, returns global min.
 * @return              The minimum locus found according to the rules above,
 *                      or `pc::MAX_IDX` if no locus exceeds the threshold.
 */
pt::idx_t find_min_locus(pt::id_t ref_id, const WalkRefIdx &wri,
                         const pvt::walk_t &w,
                         std::optional<pt::idx_t> opt_start_after) {
  bool globally = !opt_start_after.has_value();
  pt::idx_t threshold = globally ? pc::INVALID_IDX : *opt_start_after;
  pt::idx_t min_locus{pc::MAX_IDX};

  for (const auto &step : w) {
    if (auto opt_min = wri.get_min_locus(step.v_id, ref_id)) {
      if (globally || *opt_min > threshold) {
        min_locus = std::min(min_locus, *opt_min);
      }

    }
  }
return min_locus;
}

pt::idx_t comp_loop_no(pt::id_t ref_id, const WalkRefIdx &wri, const pvt::walk_t &w) {
  // the number of times the ref is seen in the walk
  // this is also the step with the max ref visits
  pt::idx_t loop_no{0};
  for (const auto &step : w) {
    loop_no = std::max(loop_no, wri.loop_count(step.v_id, ref_id));
  }

return loop_no;
}

/**
 * @brief compute min_locus and loop_no
 */
std::pair<pt::idx_t, pt::idx_t>
comp_ref_visit_bounds(pt::id_t ref_id, const WalkRefIdx &wri, const pvt::walk_t &w) {
  auto min_future = std::async(
    std::launch::async, [&]() { return find_min_locus(ref_id, wri, w, std::nullopt); });

  auto loop_future = std::async(
    std::launch::async, [&]() { return comp_loop_no(ref_id, wri, w); });

  pt::idx_t min_locus = min_future.get();
  pt::idx_t loop_no = loop_future.get();

  return {min_locus, loop_no};
}

/**
 * @brief compute itineraries for a walk
 *
 * For each ref in the walk, find the steps that are continuous with the
 * previous steps for that ref and append them to the itinerary.
 *
 * The itinerary is a collection of allele walks for each ref in the walk.
 */
void comp_itineraries(const bd::VG &g, const pvt::walk_t &w, pt::idx_t w_idx, pvt::Exp &rw) {
  // a map of ref_id to itinerary
  std::map<pt::id_t, pvt::Itn> &ref_map = rw.get_ref_itns_mut();

  WalkRefIdx wri = WalkRefIdx::from_walk(g, w);

  for (pt::id_t ref_id : wri.get_ref_ids()) { // for each ref in the walk

    auto [curr_locus, loop_count] = comp_ref_visit_bounds(ref_id, wri, w);

    for (pt::idx_t loop_no{}; loop_no < loop_count; loop_no++) {
      pvt::AW allele_walk{w_idx};

      bool is_ref_cont{false}; // is the ref continuous in the walk?

      for (pt::idx_t step_idx{}; step_idx < w.size(); ++step_idx) {

        is_ref_cont = false; // reset for each step

        const pvt::step_t &step = w[step_idx];
        pt::idx_t v_id = step.v_id;

        const VtxRefIdx &v_ref_registry = wri.get_vtx_ref_idx(v_id);

        if (v_ref_registry.has_locus(ref_id, curr_locus)) {
          {
            const bd::Vertex &v = g.get_vertex_by_id(v_id);
            std::vector<bd::RefInfo> v_ref_data = v.get_refs();

            pt::idx_t ref_data_idx = v_ref_registry.get_ref_data_idx(ref_id, curr_locus);
            auto allele_step = pvt::AS::given_ref_info(v_id, v_ref_data[ref_data_idx]);
            allele_walk.append_step(allele_step);
          }

          curr_locus += get_vtx_len(g, step);
          is_ref_cont = true; // we have a step for this ref in the walk
        }

        // we have reached a step that is not continuous with the previous
        // steps for this ref, so we should break out of the loop
        if (!is_ref_cont) {
          break;
        }
      }

      if (allele_walk.step_count() > 1) {
        pvt::Itn &itn = ref_map[ref_id];
        itn.append_at(std::move(allele_walk));
      }

      // find curr locus for next loop
      // the next loop must start at the first step in the walk
      if (loop_no + 1 < loop_count) {
        pt::idx_t nxt_curr = find_min_locus(ref_id, wri, w, curr_locus);
        if (nxt_curr == pc::INVALID_IDX) {
          break;
        }
        curr_locus = nxt_curr;
      }
    }
  }

  return;
}

} // namespace povu::genomics::allele
