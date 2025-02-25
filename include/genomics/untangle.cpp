#include "./untangle.hpp"

namespace povu::untangle {
#define MODULE "povu::untangle"

/*
 * get traversals for all refs in a single *walk*
 * assumptions:
   - a ref is contiguous in a walk
*/

std::map<pt::id_t, pvt::Itn> get_walk_refs(const bd::VG &g,
                                            const pvt::Walk &w,
                                            const std::set<pt::id_t> &ref_ids) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  std::map<pt::id_t, pvt::Itn> ref_map;

  for (auto &c_ref_id : ref_ids) {

    pvt::Itn itn {}; // associated with a single ref

    

    /* we must exclude the two boundaries otherwise any ref touching a boundary
       ends up in the walk */
    for (pt::idx_t step_idx{1}; step_idx < w.step_count() - 1; ++step_idx) {
      const pvt::Step &s = w.get_step(step_idx);
      pt::id_t v_id = s.get_v_id();
      pgt::or_e o = s.get_o();
      const bd::Vertex &v = g.get_vertex_by_id(v_id);
      const std::vector<bd::PathInfo> &v_ref_data = v.get_refs();

      pt::idx_t loop_count = 0;

      for (const bd::PathInfo &ref : v_ref_data) {
        auto [ref_id, p_o, step_idx] = ref;

        if (ref_id != c_ref_id) {
          continue;
        }

        // if (v_id > 186 && v_id < 189) {
        //   std::cerr << "v_id " << v_id << " ref_id " << ref_id << " step_idx " << step_idx << "\n";
        // }

        // TODO: do we assume that orientation is correct?
        //if (p_o != o) {
        //  std::cerr << std::format("{} WARN: walk_to_refs: mismatched orientations\n", fn_name);
        //}
        
        pvt::Step s_ = s;
        s_.set_step_idx(step_idx);

        if (itn.walk_count() <= loop_count) {
          itn.add_walk(pvt::Walk{});
        }

        itn.get_walk_mut(loop_count).append_step(std::move(s_));
        loop_count++;
      }
    }

    ref_map[c_ref_id] = std::move(itn);
  }

  return ref_map;
}

/* untangle RoVs that are flubbles linearise refs in the flubble */
pvt::RefWalks get_ref_traversals(const bd::VG &g,
                                 pvt::RoV &r,
                                 const std::set<pt::id_t> &ref_ids) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  pvt::RefWalks rw { r.get_flb() };

  //std::cerr << "\n" << "flubble " << r.get_flb().as_str() << "\n";

  /* a walk is a single traversal bounded by start to the end of an RoV */

  for (const pvt::Walk &w : r.get_walks()) {
    std::map<pt::id_t, pvt::Itn> walk_rt = get_walk_refs(g, w, ref_ids);

    // merge the ref_walks
    for (auto &[ref_id, itn] : walk_rt) {
      if (itn.walk_count() == 0) {
        continue;
      }
      rw.add_itn(ref_id, std::move(itn));
    }
  }

  rw.sort_by_step_idx();

  // print the ref walks
  // std::cerr << "\n" << "flubble " << r.get_flb().as_str() << "\n";
  // std::cerr << "tangled " << rw.is_tangled() << "\n";
  // for (const auto &[ref_id, itn] : rw.get_ref_walks()) {
  //   std::cerr << "ref " << g.get_ref_name(ref_id) << " (" << ref_id << ")\n";
  //   for (const pvt::Walk &w : itn.get_walks()) {
  //     std::cerr << w.as_str() << "\n";
  //   }
  // }
  // std::cerr << "\n";

  return rw;
}

void untangle_flb(pvt::RefWalks rt) {
  std::set<pt::id_t> ref_ids = rt.get_ref_ids();

  return;

  if (!rt.is_tangled()) {
    return;
  }

  

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

  for (pt::idx_t i {}; i < rovs.size(); ++i) {
    pvt::RoV &r = rovs[i];
    pvt::RefWalks rt = get_ref_traversals(g, r, ref_ids);
    std::cerr << "untangling " << r.get_flb().as_str() << "\n";
    untangle_flb(rt);
    std::cerr << "done untangling " << r.get_flb().as_str() << "\n";
    all_rt.push_back(rt);
  }

  // for (pvt::RoV &r : rovs) {
  //   pvt::RefWalks rt = get_ref_traversals(g, r, ref_ids);
  //   untangle_flb(rt);
  //   all_rt.push_back(rt);
  // }

  return all_rt;
}

} // namespace povu::untangle
