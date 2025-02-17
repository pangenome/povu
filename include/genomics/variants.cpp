#include "./variants.hpp"


namespace untangle {
#define MODULE "povu::untangle"

pvt::ref_walks walk_to_refs(const bd::VG &g, const pvt::Walk &w) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  pvt::ref_walks ref_map;

  for (const pvt::step &s : w.get_steps()) {
    auto [id, o] = s;
    const bd::Vertex &v = g.get_vertex_by_id(id);
    const std::vector<bd::PathInfo> &v_ref_data = v.get_refs();

    for (const bd::PathInfo &ref : v_ref_data) {
      auto [ref_id, p_o, step_idx] = ref;
      if (p_o != o) {
        std::cerr << std::format("{} WARN: walk_to_refs: mismatched orientations\n", fn_name);
      }

      ref_map[ref_id].push_back({pc::INVALID_ID, id, step_idx, o});
    }
  }

  return ref_map;
}

void merge_walk_refs(pvt::ref_walks &all_walk_refs, const pvt::ref_walks &walk_refs) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  for (const auto &[ref_id, d] : walk_refs) {
    if (all_walk_refs.find(ref_id) == all_walk_refs.end()) {
      all_walk_refs[ref_id] = d;
    } else {
      all_walk_refs[ref_id].insert(all_walk_refs[ref_id].end(), d.begin(), d.end());
    }
  }

  return;
}

void compute_loop_id(pvt::ref_walks &all_walk_refs, pvt::RoV &r) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  auto [end_id, _] = r.get_exit();
  for (auto &[ref_id, d] : all_walk_refs) {
    pt::id_t loop_id {};
    for (pvt::ref_step_t &step : d) {
     step.set_loop_id(loop_id);
     if (step.get_v_id() == end_id) {
       loop_id++;
     }
    }
  }

  return;
}

/* untangle RoVs that are flubbles */
void untangle_flb(const bd::VG &g, pvt::RoV &r) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  pvt::ref_walks all_walk_refs;

  for (const pvt::Walk &w : r.get_walks()) {
    pvt::ref_walks walk_refs = walk_to_refs(g, w);
    merge_walk_refs(all_walk_refs, walk_refs);
  }

  // sort by step_idx
  for (auto &[ref_id, d] : all_walk_refs) {
    std::sort(d.begin(), d.end(), [](const pvt::ref_step_t &a, const pvt::ref_step_t &b) {
      return a.step_idx < b.step_idx;
    });
  }

  compute_loop_id(all_walk_refs, r);

  return;
}

} // namespace untangle

namespace povu::variants {
#define MODULE "povu::variants"

std::vector<pvt::RoV> par_populate_walks(const bd::VG &g, std::vector<pvt::RoV> rs) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  std::size_t cfl_count = rs.size();

  uint8_t num_threads = std::thread::hardware_concurrency();
  if (num_threads == 0) {
    num_threads = 1;
  }
  // avoid creating more threads than flubbles
  if (static_cast<std::size_t>(num_threads) > cfl_count) {
    num_threads = cfl_count;
  }

  std::vector<std::thread> threads;
  threads.reserve(num_threads);

  //std::size_t fl_count = can_fls.size();
  std::size_t chunk_size = cfl_count / num_threads;
  std::size_t remainder = cfl_count % num_threads;
  std::size_t start { 0 };

  auto worker = [&](std::size_t start, std::size_t end) {
    for (std::size_t i = start; i < end; ++i) {
      bd::populate_walks(g, rs[i], MAX_FLUBBLE_STEPS);
    }
  };

  for (uint8_t i = 0; i < num_threads; ++i) {
    std::size_t end = start + chunk_size + (i < remainder ? 1 : 0);
    threads.push_back(std::thread(worker, start, end));
    start = end;
  }

  for (auto &t : threads) {
    t.join();
  }

  return rs;
}

 /* filter out RoVs whose walk count is less than 2 */
inline void filter_invalid(std::vector<pvt::RoV> &r) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};
  for (auto it = r.begin(); it != r.end(); ++it) {
    if (it->walk_count() < 2) {
      std::cerr << std::format("{} WARN: flubble {} has {} paths\n",
                               fn_name, it->as_str(), it->walk_count());
      it = r.erase(it);
    }
  }
};

/* initialize RoVs from flubbles */
inline std::vector<pvt::RoV> init_rovs(const std::vector<pgt::flubble_t> &fls) {
  std::vector<pvt::RoV> rs;
  rs.reserve(fls.size());
  for (auto &fl : fls) {
    rs.push_back(pvt::RoV{fl});
  }

  return rs;
}


void call_variants(const std::vector<pgt::flubble_t>& canonical_flubbles,
                   const bd::VG& g,
                   const core::config& app_config) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  std::vector<pvt::RoV> rs = init_rovs(canonical_flubbles);
  std::vector<pvt::RoV> rovs = par_populate_walks(g, rs); // par for parallel
  filter_invalid(rovs);

  for (auto &rov : rovs) {
    untangle::untangle_flb(g, rov);
  }
}


} // namespace povu::variants
