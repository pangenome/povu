#include "./untangle.hpp"

namespace povu::untangle {
#define MODULE "povu::untangle"

/* get traversals for all refs in a single *walk* */
std::map<pt::id_t, pvt::Walk> get_walk_refs(const bd::VG &g,
                                            const pvt::Walk &w,
                                            const std::set<pt::id_t> &ref_ids) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  std::map<pt::id_t, pvt::Walk> ref_map;

  /* we must exclude the two boundaries otherwise any ref touching a bounary ends up in the walk */

  for (pt::idx_t step_idx {1}; step_idx < w.step_count() - 1; ++step_idx) {
    const pvt::Step &s = w.get_step(step_idx);
    pt::id_t v_id = s.get_v_id();
    pgt::or_e o = s.get_o();
    const bd::Vertex &v = g.get_vertex_by_id(v_id);
    const std::vector<bd::PathInfo> &v_ref_data = v.get_refs();

    for (const bd::PathInfo &ref : v_ref_data) {
      auto [ref_id, p_o, step_idx] = ref;
      if (!ref_ids.contains(ref_id)) {
        continue;
      }

      if (p_o != o) {
        std::cerr << std::format(
            "{} WARN: walk_to_refs: mismatched orientations\n", fn_name);
      }

      pvt::Step s_ = s;
      s_.set_step_idx(step_idx);
      ref_map[ref_id].append_step(std::move(s_));
    }
  }

  //for (const pvt::Step &s : w.get_steps()) {

  //}

    return ref_map;
}

/* untangle RoVs that are flubbles linearise refs in the flubble */
pvt::RefWalks get_ref_traversals(const bd::VG &g,
                                 pvt::RoV &r,
                                 const std::set<pt::id_t> &ref_ids) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  pvt::RefWalks rw { r.get_flb() };

  /* a walk is a single traversal from the start to the end of an RoV */

  for (const pvt::Walk &w : r.get_walks()) {
    std::map<pt::id_t, pvt::Walk> walk_rt = get_walk_refs(g, w, ref_ids);

    std::cerr << fn_name << " " << w.as_str() << std::endl;

    
    // merge the ref_walks
    for (auto &[ref_id, w] : walk_rt) {
      std::cerr << "ref " << g.get_ref_name(ref_id) << " (" << ref_id << ")\n";
      rw.add_walk(ref_id, std::move(w));
    }
  }

  rw.sort_by_step_idx();

  return rw;
}

void untangle_flb(pvt::RefWalks rt) {
  std::set<pt::id_t> ref_ids = rt.get_ref_ids();

  // all vs all compare
 for (pt::id_t ref_id1 : ref_ids) {
   const pvt::Itn &rw1 = rt.get_itn(ref_id1);
   for (pt::id_t ref_id2 : ref_ids) {
     if (ref_id1 == ref_id2 || rt.has_aln(ref_id1, ref_id2)) {
       continue;
     }
     const pvt::Itn &rw2 = rt.get_itn(ref_id2);

     // compare at the at level
     std::string et = pa::align(rw1, rw2, pvt::aln_level_e::at);

     rt.add_aln(ref_id1, ref_id2, std::move(et));
   }
 }

 return;
}

std::vector<pvt::RefWalks> untangle_flb_rovs(const bd::VG &g,
                                             std::vector<pvt::RoV> &rovs,
                                             const std::set<pt::id_t> &ref_ids) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  std::vector<pvt::RefWalks> all_rt;
  all_rt.reserve(rovs.size());
  for (pvt::RoV &r : rovs) {
    pvt::RefWalks rt = get_ref_traversals(g, r, ref_ids);
    untangle_flb(rt);
    all_rt.push_back(rt);
  }

  return all_rt;
}

} // namespace povu::untangle
