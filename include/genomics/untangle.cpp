#include "./untangle.hpp"
#include <sys/types.h>
#include <unordered_set>


namespace povu::untangle {
#define MODULE "povu::untangle"

// typedef std::map<pt::id_t, pvt::Itn> RefItns;
//typedef std::pair<pt::id_t, pt::id_t> up_t;

/**
 * get traversals for all refs in a single *walk*
 * assumptions:
 * - a ref is contiguous in a walk
 * - orientation between vertex and ref is correct
*/
std::map<pt::id_t, pvt::Itn> get_itns(const bd::VG &g, const pvt::Walk &w) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  std::map<pt::id_t, pvt::Itn> ref_map;

  /* we must exclude the two boundaries otherwise any ref touching a boundary
     ends up in the walk */
  for (pt::idx_t step_idx {}; step_idx < w.step_count(); ++step_idx) {

    pt::id_t v_id = w.get_step(step_idx).get_v_id();
    const bd::Vertex &v = g.get_vertex_by_id(v_id);
    std::vector<bd::PathInfo> v_ref_data = v.get_refs();

    std::map<pt::id_t, pt::idx_t> loop_count;

    for (const bd::PathInfo &ref : v_ref_data) {
      auto [ref_id, p_o, locus] = ref;

      // Using operator[] initializes lc to 0 if not present.
      pt::idx_t &lc = loop_count[ref_id];

      // Try to insert a new Itn if one does not exist for ref_id.
      auto [it, inserted] = ref_map.try_emplace(ref_id, pvt::Itn{});
      pvt::Itn &itn = it->second;

      if (itn.at_count() <= lc) {
        itn.add_at(pvt::AT{});
      }

      if (itn.get_at(lc).step_count() > 0) {
        const pvt::Step prev_step = itn.get_at(lc).get_steps().back();
        pt::idx_t prev_locus = prev_step.get_step_idx();
        pt::id_t prev_v_id = prev_step.get_v_id();
        const bd::Vertex &prev_v = g.get_vertex_by_id(prev_v_id);
        pt::idx_t vtx_len = prev_v.get_label().length();

        if (locus !=  vtx_len + prev_locus) {
          continue;
        }
      }


      itn.append_step(lc, pvt::Step{v_id, locus, p_o});
      lc++;
    }
  }

  // how to avoid this
  // because we only pass those with more than 2 steps
  // delete those with a single step
  for (auto it = ref_map.begin(); it != ref_map.end();) {
    if (it->second.step_count() < 3) {
      it = ref_map.erase(it);
    } else {
      ++it;
    }
  }


  return ref_map;
}

/**
  assumption:
   - a deletion is in both s and t
     * even if orientation matches RoV (s or t) and is only in s or t it is not
       a deletion
   - they have the same loop count in both s and t
*/
std::map<pt::id_t, pvt::Itn> get_del_itns(const bd::VG &g, pvt::RoV &r) {

  auto [s_v_id, _] = r.get_entry();
  const bd::Vertex &s_vtx = g.get_vertex_by_id(s_v_id);
  std::vector<bd::PathInfo> s_ref_data = s_vtx.get_refs();
  pt::idx_t s_vtx_len = s_vtx.get_label().length();

  auto [t_v_id, __] = r.get_exit();
  const bd::Vertex &t_vtx = g.get_vertex_by_id(t_v_id);
  std::vector<bd::PathInfo> t_ref_data = t_vtx.get_refs();

  std::map<pt::id_t, pvt::Itn> ref_map;

  for (auto [ref_id, o, locus] : s_ref_data) {
    for (auto [ref_id2, o2, locus2] : t_ref_data) {
      if (ref_id == ref_id2 && locus + s_vtx_len == locus2) {

        pvt::AT at;
        at.append_step(pvt::Step{s_v_id, locus, o});
        at.append_step(pvt::Step{t_v_id, locus2, o2});
        ref_map[ref_id].add_at(std::move(at));
      }
    }
  }

  return ref_map;
}

/* untangle RoVs that are flubbles linearise refs in the flubble */
pvt::RefWalks get_ref_traversals(const bd::VG &g, pvt::RoV &r) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  pvt::RefWalks rw { r.get_flb() };

  // merge the allele traversals
  auto merge_ats = [&](std::map<pt::id_t, pvt::Itn> walk_rt) {
    for (auto &[ref_id, itn] : walk_rt) {
      rw.add_itn(ref_id, std::move(itn));
    }
  };

  // a walk is a single traversal bounded by start to the end of an RoV
  std::map<pt::id_t, pvt::Itn> walk_rt;
  for (const pvt::Walk &w : r.get_walks()) {

    // should this be here?
    if (w.step_count() < 3) {
      continue;
    }

    walk_rt = get_itns(g, w);
    merge_ats(walk_rt);
  }

  // merge the del allele traversals
  walk_rt = get_del_itns(g, r);
  merge_ats(walk_rt);

  rw.sort_by_step_idx();


  // print the ref walks
  std::cerr << "flubble " << r.get_flb().as_str() << "\n";
  std::cerr << "is tangled " << rw.is_tangled() << "\n";
  for (auto &[ref_id, itn] : rw.get_ref_walks()) {
    std::cerr << "ref " << g.get_ref_name(ref_id) << "\n";
    for (auto &at : itn.get_ats()) {
      std::cerr << at.as_str() << "\n";
    }
  }

  return rw;
}

inline std::vector<pt::up_t<pt::id_t>> compute_pairs(pvt::RefWalks rt) {
  std::set<pt::id_t> ref_ids = rt.get_ref_ids();

  std::set<pt::up_t<pt::id_t>> done;
  std::vector<pt::up_t<pt::id_t>> aln_pairs;

  // all vs all compare
  for (pt::id_t ref_id1 : ref_ids) {
    for (pt::id_t ref_id2 : ref_ids) {

      if (ref_id1 == ref_id2) {
        continue;
      }

      pt::up_t<pt::id_t> p{ref_id1, ref_id2};

      if (done.find(p) != done.end()) {
        continue;
      }

      done.insert(p);
      aln_pairs.push_back(p);
    }
  }

  return aln_pairs;
}

void untangle_flb(pvt::RefWalks rt) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  std::vector<pt::up_t<pt::id_t>> aln_pairs = compute_pairs(rt);

  for (auto [ref_id1, ref_id2] : aln_pairs) {

    const pvt::Itn &itn1 = rt.get_itn(ref_id1);
    const pvt::Itn &itn2 = rt.get_itn(ref_id2);

    // std::cerr << "ref pairs " << g.get_ref_name(ref_id1) << " vs " << g.get_ref_name(ref_id2) << "\n";
    // std::cerr << itn1.at_count() << " " << itn2.at_count() << "\n";
    // for (auto x : itn1.get_ats()) {
    //   std::cerr << x.as_str() << ", ";
    // }
    // std::cerr << " vs  ";
    // for (auto x : itn2.get_ats()) {
    //   std::cerr << x.as_str() << ", ";
    // }
    // std::cerr << "\n";


    std::string et = pa::align(itn1, itn2, pvt::aln_level_e::at);

    //std::cerr << et << "\n";

    rt.add_aln(ref_id1, ref_id2, std::move(et));
  }
}

std::vector<pvt::RefWalks> untangle_flb_rovs(const bd::VG &g, std::vector<pvt::RoV> &rovs) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  std::vector<pvt::RefWalks> all_rt;
  all_rt.reserve(rovs.size());

  for (pt::idx_t i {}; i < rovs.size(); ++i) {
    pvt::RoV &r = rovs[i];
    pvt::RefWalks rt = get_ref_traversals(g, r);

    if (rt.is_tangled()) {
      untangle_flb(rt);
    }

    all_rt.push_back(rt);
  }

  return all_rt;
}

} // namespace povu::untangle
