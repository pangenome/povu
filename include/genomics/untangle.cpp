#include "./untangle.hpp"
#include <vector>


namespace povu::untangle {
#define MODULE "povu::untangle"

pvt::ref_walks get_ref_walks(const bd::VG &g, const pvt::Walk &w) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  pvt::ref_walks ref_map;

  for (const pvt::Step &s : w.get_steps()) {
    //auto [id, o] = s;
    pt::id_t v_id = s.get_v_id();
    pgt::or_e o = s.get_o();
    const bd::Vertex &v = g.get_vertex_by_id(v_id);
    const std::vector<bd::PathInfo> &v_ref_data = v.get_refs();

    for (const bd::PathInfo &ref : v_ref_data) {
      auto [ref_id, p_o, step_idx] = ref;
      if (p_o != o) {
        std::cerr << std::format("{} WARN: walk_to_refs: mismatched orientations\n", fn_name);
      }

      pvt::Step s_ = s;
      s_.set_step_idx(step_idx);
      ref_map[ref_id].push_back(s_);
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

void set_loop_id(pvt::ref_walks &all_walk_refs, pvt::RoV &r) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  auto [end_id, _] = r.get_exit();
  for (auto &[ref_id, rw] : all_walk_refs) {
    pt::id_t loop_id {};
    for (pvt::Step &s : rw) {
     s.set_loop_id(loop_id);
     if (s.get_v_id() == end_id) {
       loop_id++;
     }
    }
  }

  return;
}

/* untangle RoVs that are flubbles linearise refs in the flubble */
pvt::ref_walks linearise_flb_refs(const bd::VG &g, pvt::RoV &r) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  pvt::ref_walks all_walk_refs;

  for (const pvt::Walk &w : r.get_walks()) {
    pvt::ref_walks walk_refs = get_ref_walks(g, w);
    merge_walk_refs(all_walk_refs, walk_refs);
  }

  // sort by step_idx within each ref
  for (auto &[ref_id, d] : all_walk_refs) {
    std::sort(d.begin(), d.end(), [](const pvt::Step &a, const pvt::Step &b) {
      return a.get_step_idx() < b.get_step_idx();
    });
  }

  set_loop_id(all_walk_refs, r);

  return all_walk_refs;
}

void compare_refs(const pvt::ref_walk &iw, const pvt::ref_walk &jw) {

  // align step/node level
  pa::align_steps(iw, jw);

  // align RoV level
  pa::align_rovs(iw, jw);

}

void untangle_flb(pvt::ref_walks rws) {
  std::size_t ref_count = rws.size();

  for (std::size_t i {}; i < ref_count; i++) {
    for (std::size_t j {}; j < ref_count; j++) {
      if (i == j) {
        continue;
      }

      pvt::ref_walk iw = rws[i];
      pvt::ref_walk jw = rws[j];

      compare_refs(iw, jw);

    }
  }

}

void untangle_flb_rovs(const bd::VG &g, std::vector<pvt::RoV> rovs) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  std::vector<pvt::ref_walks> linearised_refs;
  linearised_refs.reserve(rovs.size());
  for (pvt::RoV &r : rovs) {
    pvt::ref_walks rw = linearise_flb_refs(g, r);
    linearised_refs.push_back(rw);
  }

  /*  */
  for (auto rws : linearised_refs) {
    untangle_flb(rws);
  }

  return;
}

} // namespace povu::untangle
