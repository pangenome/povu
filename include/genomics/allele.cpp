#include "./allele.hpp"
#include <cstdlib>
#include <map>
#include <optional>
#include <sys/types.h>
#include <thread>
#include <unordered_set>
#include <vector>

namespace povu::genomics::allele {


pt::idx_t get_vtx_len(const bd::VG &g, const pgt::step_t &s) {
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
pt::idx_t find_min_locus(const bd::VG &g, pt::id_t ref_id, const pgt::walk_t &w,
                         std::optional<pt::idx_t> opt_start_after) {
  pt::idx_t min_locus{pc::MAX_IDX};

  for (const auto &[v_id, _] : w) {
    const pgr::VtxRefData &v_ref_data = g.get_vertex_by_id(v_id).get_refs();
    //const bd::VtxRefIdx &vr_idx = g.get_vtx_ref_idx(v_id);
    if (auto opt_min = v_ref_data.get_min_locus(ref_id, opt_start_after)) {
      min_locus = std::min(min_locus, *opt_min);
    }
  }

  return min_locus;
}


pt::idx_t comp_loop_no(const bd::VG &g, pt::id_t ref_id, const pgt::walk_t &w) {
  // the number of times the ref is seen in the walk
  // this is also the step with the max ref visits
  pt::idx_t max_loop_no{0};
  for (const auto &[v_id, _] : w) {
    const pgr::VtxRefData &v_ref_data = g.get_vertex_by_id(v_id).get_refs();
    pt::idx_t times = v_ref_data.loop_count(ref_id);
    max_loop_no = std::max(max_loop_no, v_ref_data.loop_count(ref_id));
  }

return max_loop_no;
}

/**
 * @brief compute min_locus and loop_no
 */
std::pair<pt::idx_t, pt::idx_t> comp_ref_visit_bounds(const bd::VG &g,
                                                      pt::id_t ref_id,
                                                      const pgt::walk_t &w) {
  auto min_future = std::async(
    std::launch::async, [&]() { return find_min_locus(g, ref_id, w, std::nullopt); });

  auto loop_future = std::async(
    std::launch::async, [&]() { return comp_loop_no(g, ref_id, w); });

  pt::idx_t min_locus = min_future.get();
  pt::idx_t loop_no = loop_future.get();

  return {min_locus, loop_no};
}

inline std::optional<AW> do_walk(const bd::VG &g, pt::id_t ref_id, pt::idx_t w_idx,
                          const pgt::walk_t &w, pt::idx_t start_locus) {
  pt::idx_t curr_locus = start_locus;
  AW allele_walk{w_idx};

  for (pt::idx_t step_idx{}; step_idx < w.size(); ++step_idx) {

    const pgt::step_t &step = w[step_idx];
    auto [v_id, o] = step;

    const bd::Vertex v = g.get_vertex_by_id(v_id);
    const pgr::VtxRefData &ref_data = v.get_refs();

    pgr::ref_step_t ref_addr = {o, curr_locus};
    if (!ref_data.has_address(ref_id, ref_addr)) {
      // we have reached a step that is not continuous with the previous
      // steps for this ref, so we should break out of the loop
      break;
    }

    auto allele_step = AS::given_ref_info(v_id, ref_addr);
    allele_walk.append_step(std::move(allele_step));

    curr_locus += v.get_label().length();
  }

  if (allele_walk.step_count() <= 1) {
    return std::nullopt;
  }

  return allele_walk;
};

std::vector<AW> handle_walks(const bd::VG &g, pt::id_t ref_id, pt::idx_t w_idx,
                             const pgt::walk_t &w, pt::idx_t start_locus,
                             pt::idx_t loop_count) {
  std::vector<AW> aws;
  aws.reserve(loop_count);
  if (std::optional<AW> opt_aw = do_walk(g, ref_id, w_idx, w, start_locus)) {
    aws.push_back(std::move(*opt_aw));
  }

  if (loop_count == 1) {
    return aws;
  }

  auto [first_v_id, __] = w.front(); // first step in the walk
  const pgr::VtxRefData &v_ref_data = g.get_vertex_by_id(first_v_id).get_refs();
  const std::set<pgr::ref_step_t> &starts = v_ref_data.get_ref_data(ref_id);

  for (auto [_, loop_start_locus] : starts) {
    if (loop_start_locus != start_locus) {
      if (std::optional<AW> opt_aw = do_walk(g, ref_id, w_idx, w, loop_start_locus)) {
        aws.push_back(std::move(*opt_aw));
      }
    }
  }

  return aws;
}

std::vector<std::pair<pt::idx_t, id_t>>
pre_comp_tasks(const bd::VG &g, const std::vector<pgt::walk_t> &walks ) {
  std::vector<std::pair<pt::idx_t, id_t>> tasks;
  for (pt::idx_t w_idx = 0; w_idx < walks.size(); ++w_idx) {
    const pgt::walk_t &w = walks[w_idx];
    std::set<pt::id_t> ref_ids;

    for (const pgt::step_t &s : w) {
      auto v_id = s.v_id;
      const std::set<pt::idx_t> &step_ref_ids = g.get_vertex_by_id(v_id).get_refs().get_ref_ids();
      ref_ids.insert(step_ref_ids.begin(), step_ref_ids.end());
    }

    for (pt::id_t ref_id : ref_ids) {
      tasks.push_back({w_idx, ref_id});
    }
  }

  return tasks;
}

/**
 * @brief compute itineraries for a walk
 *
 * For each ref in the walk, find the steps that are continuous with the
 * previous steps for that ref and append them to the itinerary.
 *
 * The itinerary is a collection of allele walks for each ref in the walk.
 */
void comp_itineraries_async(const bd::VG &g,
                            const std::vector<pgt::walk_t> &walks,
                            std::map<pt::id_t, Itn> &ref_map,
                            povu::thread::thread_pool &pool) {

  std::vector<std::pair<pt::idx_t, pt::id_t>> tasks = pre_comp_tasks(g, walks);
  const std::size_t M = tasks.size();

  std::vector<std::future<std::vector<AW>>> futs;
  futs.reserve(M);

  for (auto [w_idx, ref_id] : tasks) {
    futs.emplace_back(pool.submit([&, w_idx, ref_id]() -> std::vector<AW> {
      const auto &w = walks[w_idx];
      auto [start_locus, loop_count] = comp_ref_visit_bounds(g, ref_id, w);
      return handle_walks(g, ref_id, w_idx, w, start_locus, loop_count);
    }));
  }

  // Merge results on this thread (no locks needed)
  for (std::size_t i = 0; i < M; ++i) {
    auto aws = futs[i].get(); // rethrows from worker
    if (aws.empty())
      continue;
    auto ref_id = tasks[i].second;
    Itn &itn = ref_map[ref_id];
    for (AW &aw : aws)
      itn.append_at(std::move(aw));
  }

  return;
}

void comp_itineraries_serial(const bd::VG &g, const std::vector<pgt::walk_t> &walks,
                      std::map<pt::id_t, Itn> &ref_map) {

  std::vector<std::pair<pt::idx_t, id_t>> tasks = pre_comp_tasks(g, walks);

  for (auto [w_idx, ref_id] : tasks) {
    const auto &w = walks[w_idx];
    auto [start_locus, loop_count] = comp_ref_visit_bounds(g, ref_id, w);
    std::vector<AW> aws = handle_walks(g, ref_id, w_idx, w, start_locus, loop_count);

    if (!aws.empty()) {
      Itn &itn = ref_map[ref_id];
      for (AW &aw : aws) {
        itn.append_at(std::move(aw));
      }
    }
  }

  return;
}

void comp_itineraries(const bd::VG &g, const std::vector<pgt::walk_t> &walks,
                      std::map<pt::id_t, Itn> &ref_map,
                      povu::thread::thread_pool &pool) {

  if (walks.size() > 10) {
    comp_itineraries_async(g, walks, ref_map, pool);
  }
  else {
    comp_itineraries_serial(g, walks, ref_map);
  }

  return;
}

} // namespace povu::genomics::allele
