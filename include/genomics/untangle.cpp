#include "./untangle.hpp"
#include <string>
#include <vector>


namespace povu::untangle {
#define MODULE "povu::untangle"

/* get traversals for all refs in a single *walk* */
std::map<pt::id_t, pvt::Walk> get_walk_refs(const bd::VG &g, const pvt::Walk &w) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  std::map<pt::id_t, pvt::Walk> ref_map;

  for (const pvt::Step &s : w.get_steps()) {
    pt::id_t v_id = s.get_v_id();
    pgt::or_e o = s.get_o();
    const bd::Vertex &v = g.get_vertex_by_id(v_id);
    const std::vector<bd::PathInfo> &v_ref_data = v.get_refs();

    for (const bd::PathInfo &ref : v_ref_data) {
      auto [ref_id, p_o, step_idx] = ref;
      if (p_o != o) {
        std::cerr << std::format(
            "{} WARN: walk_to_refs: mismatched orientations\n", fn_name);
      }

      pvt::Step s_ = s;
      s_.set_step_idx(step_idx);
      ref_map[ref_id].append_step(std::move(s_));
    }
  }

  return ref_map;
}

/* untangle RoVs that are flubbles linearise refs in the flubble */
pvt::RefWalks get_ref_traversals(const bd::VG &g, pvt::RoV &r) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  pvt::RefWalks rw;

  /* a walk is a single traversal from the start to the end of an RoV */

  for (const pvt::Walk &w : r.get_walks()) {
    std::map<pt::id_t, pvt::Walk> walk_rt = get_walk_refs(g, w);
    // merge the ref_walks
    for (auto &[ref_id, w] : walk_rt) {
      rw.add_walk(ref_id, std::move(w));
    }
  }

  rw.sort_by_step_idx();
  return rw;
}


void compare_traversals(const pvt::Itn &rw1, const pvt::Itn &rw2) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  // compare at the itn level
  std::string et = pa::align(rw1, rw2, pvt::aln_level_e::rov);

  bool aln_step = false;
  for (auto c : et) {
    if (c != 'M') {
      aln_step= true;
      break;
    }
  }

  std::string et_step;
  if (aln_step) {
    et_step = pa::align(rw1, rw2, pvt::aln_level_e::step);
  }

  return;
}

void untangle_flb(pvt::RefWalks rt) {

  // up = unordered pair
  typedef std::pair<pt::id_t, pt::id_t> up_t;
  std::set<up_t> compared;

  auto to_up = [](pt::id_t a, pt::id_t b) -> up_t {
    return {std::min(a, b), std::max(a, b)};
  };

  auto add_to_compared = [&](pt::id_t ref_id1, pt::id_t ref_id2) {
    compared.insert(to_up(ref_id1, ref_id2));
  };

  auto is_compared = [&](pt::id_t ref_id1, pt::id_t ref_id2) -> bool {
    return compared.contains(to_up(ref_id1, ref_id2));
  };

  std::set<pt::id_t> ref_ids = rt.get_ref_ids();

  // all vs all compare
 for (pt::id_t ref_id1 : ref_ids) {
   const pvt::Itn &rw1 = rt.get_itn(ref_id1);
   for (pt::id_t ref_id2 : ref_ids) {
     if (ref_id1 == ref_id2 || is_compared(ref_id1, ref_id2)) {
       continue;
     }
     const pvt::Itn &rw2 = rt.get_itn(ref_id2);
     compare_traversals(rw1, rw2);
     add_to_compared(ref_id1, ref_id2);
   }
 }

 return;
}

void untangle_flb_rovs(const bd::VG &g, std::vector<pvt::RoV> rovs) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  for (pvt::RoV &r : rovs) {
    pvt::RefWalks rt = get_ref_traversals(g, r);
    untangle_flb(rt);
  }

  return;
}

} // namespace povu::untangle
