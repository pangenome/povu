#include "./untangle.hpp"
#include <format>
#include <utility>
#include <vector>


namespace povu::untangle {
#define MODULE "povu::untangle"


std::map<pt::id_t, pt::idx_t> count_max_ats(const bd::VG &g, const pvt::Walk &w) {

  std::map<pt::id_t, pt::idx_t> loop_count;

  /*
   * we must exclude the two boundaries otherwise any ref touching a boundary
   * ends up in the walk
   */
  for (pt::idx_t step_idx{}; step_idx < w.step_count(); ++step_idx) {

    pt::id_t v_id = w.get_step(step_idx).get_v_id();
    const bd::Vertex &v = g.get_vertex_by_id(v_id);
    std::vector<bd::PathInfo> v_ref_data = v.get_refs();

    // ref id to loop count
    std::map<pt::id_t, pt::idx_t> local_loop_count;

    // count the number of times a ref is seen in a step
    for (const bd::PathInfo &ref : v_ref_data) {
      auto [ref_id, p_o, locus] = ref;

      local_loop_count[ref_id]++;
    }

    // set the max loop count for each ref
    for (auto &[ref_id, lc] : local_loop_count) {
      if (loop_count.find(ref_id) == loop_count.end()) {
        loop_count[ref_id] = lc;
      } else {
        if (lc > loop_count[ref_id]) {
          loop_count[ref_id] = lc;
        }
      }
    }
  }

  return loop_count;
}

/**
 * get traversals for all refs in a single *walk*
 * assumptions:
 *  - a ref is contiguous in a walk
 *  - orientation between vertex and ref is correct
 *  - a ref will not show up in only one vertex in the walk
 * - a deletion is in both s and t
     * even if orientation matches RoV (s or t) and is only in s or t it is not
       a deletion
   - they have the same loop count in both s and t
 * @param g the variation graph
 * @param w a single traversal bounded by start to the end of an RoV
 * @param rws the ref walks
 */
void get_itns(const bd::VG &g, pt::idx_t walk_idx, const pvt::Walk &w, pvt::RefWalks &rws) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  pgt::flubble_t fl = rws.get_flb();
  auto [s_v_id, _] = fl.start_;
  auto [e_v_id, __] = fl.end_;


  std::map<pt::id_t, pt::idx_t> loop_count = count_max_ats(g, w);

  std::map<pt::id_t, pvt::Itn> ref_map;
  /* init ref map */
  for (auto &[ref_id, lc] : loop_count) {
    ref_map[ref_id] = pvt::Itn{};
    for (pt::idx_t i{}; i < lc; ++i) {
      ref_map[ref_id].append_at(pvt::AT{});
    }
  }

  auto set_deletions = [&](std::map<pt::id_t, pvt::Itn> &ref_map) {
    for (auto &[ref_id, itn] : ref_map) {
      // mark deletions as deletions
      for (auto &at : itn.get_ats_mut()) {
        if (at.step_count() == 2 &&
            at.get_steps()[0].get_v_id() == s_v_id &&
            at.get_steps()[1].get_v_id() == e_v_id) {
          at.set_is_del(true);
        }
      }
    }
  };

  // merge the allele traversals
  auto merge_ats = [&](std::map<pt::id_t, pvt::Itn> &walk_rt) {
    for (auto &[ref_id, itn] : walk_rt) {
      rws.add_itn(ref_id, std::move(itn));
    }
  };

  auto is_extendable = [&](const pvt::AT &curr_at, pt::idx_t step_idx, pt::idx_t locus) -> bool {
    const std::vector<pvt::Step> &steps = curr_at.get_steps();

    const pvt::Step &s = steps.back();
    pt::idx_t prev_locus = s.get_step_idx();
    pt::id_t prev_v_id = s.get_v_id();
    const bd::Vertex &prev_v = g.get_vertex_by_id(prev_v_id);
    pt::idx_t vtx_len = prev_v.get_label().length();

    pt::id_t prev_step_v_id = w.get_step(step_idx - 1).get_v_id();

    return (prev_v_id == prev_step_v_id && prev_locus + vtx_len == locus);
  };

  /*
   * we must exclude the two boundaries otherwise any ref touching a
   * boundary ends up in the walk
   */
  for (pt::idx_t step_idx{}; step_idx < w.step_count(); ++step_idx) {

    pt::id_t v_id = w.get_step(step_idx).get_v_id();
    const bd::Vertex &v = g.get_vertex_by_id(v_id);
    std::vector<bd::PathInfo> v_ref_data = v.get_refs();

    for (const bd::PathInfo &ref : v_ref_data) {
      auto [ref_id, p_o, locus] = ref;

      // std::cerr << "ref " << g.get_ref_name(ref_id) << "\n";

      pvt::Itn &curr_itn = ref_map[ref_id];

      // if the start of a new AT then find the first empty itn and add
      // the step it is a start if:
      //  - the step idx is 0
      //  - step_idx > 0 but cannot extend

      /* attempt to extend an at */
      std::vector<pvt::AT> &ats = curr_itn.get_ats_mut();
      for (pt::idx_t at_idx{}; at_idx < ats.size(); ++at_idx) {
        pvt::AT &curr_at = ats[at_idx];

        if (curr_at.step_count() == 0 && step_idx == 0) {
          curr_at.append_step(pvt::Step{v_id, locus, p_o});
          break;
        }

        if (step_idx == 0) {
          continue;
        }

        if (is_extendable(curr_at, step_idx, locus)) {
          curr_at.append_step(pvt::Step{v_id, locus, p_o});
          break;
        }
      }
    }
  }

  set_deletions(ref_map);
  merge_ats(ref_map);
}

pvt::Itn remove_invalid_ats(pvt::Itn &itn) {
  //pt::idx_t curr_locus { pc::INVALID_IDX };

  // locus and index at idx of the longest AT for that locus
  std::map<pt::id_t, std::pair<pt::idx_t, pt::idx_t>> longest;

  //std::vector<pt::idx_t> to_del;

  auto foo = [&](pt::idx_t at_idx) {
    const pvt::AT &at = itn.get_at(at_idx);
    pt::idx_t sc = at.step_count();
    pt::idx_t lc = at.get_steps().front().get_step_idx();

    return std::make_pair(sc, lc);
  };

  for (pt::idx_t at_idx{}; at_idx < itn.at_count(); at_idx++) {
    auto [sc, lc] = foo(at_idx);

    if (longest.contains(lc) && sc > longest[lc].first) {
      longest[lc] = std::make_pair(sc, at_idx);
    }
    else if (!longest.contains(lc)) {
      longest[lc] = std::make_pair(sc, at_idx);
    }
  }

  // create a new itn with the longest AT per locus
  pvt::Itn itn_ {};
  for (auto [_, v] : longest) {
    auto [__, at_idx] = v;
    if (at_idx >= itn.at_count()) {
      auto s = std::format("at index {} out of range {}", at_idx, itn.at_count());
      throw std::out_of_range(s);
    }
    pvt::AT at = itn.get_at(at_idx);
    itn_.append_at(std::move(at));
  }

  return itn_;


}

/* untangle RoVs that are flubbles linearise refs in the flubble */
pvt::RefWalks get_ref_traversals(const bd::VG &g, pvt::RoV &r) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  pvt::RefWalks rw { r.get_flb() };

  /* a walk is a single traversal bounded by start to the end of an RoV */
  for (pt::idx_t walk_idx{}; walk_idx < r.walk_count(); ++walk_idx) {
    pvt::Walk const &w = r.get_walks()[walk_idx];
    get_itns(g, walk_idx, w, rw);
  }


  rw.sort_by_step_idx();

  for (auto &[r, itn] : rw.get_ref_itns_mut()) {
    pvt::Itn n = remove_invalid_ats(itn);
    rw.replace_itn(r, std::move(n));
  }

  for (auto &[_, itn] : rw.get_ref_itns()) {
    if (itn.at_count() > 1) {
      rw.set_tangled(true);
      break;
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

/**
 * align the traversals of two refs
 */
void untangle_flb(const bd::VG &g, pvt::RefWalks &rt) {
  std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  std::vector<pt::up_t<pt::id_t>> aln_pairs = compute_pairs(rt);

  for (auto [ref_id1, ref_id2] : aln_pairs) {

    const pvt::Itn &itn1 = rt.get_itn(ref_id1);
    const pvt::Itn &itn2 = rt.get_itn(ref_id2);

    std::string et = pa::align(itn1, itn2, pvt::aln_level_e::at);

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
      untangle_flb(g, rt);
    }

    all_rt.push_back(rt);
  }

  return all_rt;
}

} // namespace povu::untangle
