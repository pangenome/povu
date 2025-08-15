#include "./utils.hpp"

namespace povu::genomics::utils {


// should this be genotyping?
namespace variants {


// the ref visits of a single vertex. Unlike the ref visits in a walk, this
// contains the steps that the ref takes in the vertex, sorted by step index
typedef std::vector<pvt::AS> VtxRefVisits;


class VtxRefMeta {
  // a map of ref id to steps it has in the vertex
  std::map<pt::id_t, VtxRefVisits> data_;

public:
  // --------------
  // constructor(s)
  // --------------

  VtxRefMeta(){}

  // ---------
  // getter(s)
  // ---------

  std::set<pt::id_t> get_refs() const {
    std::set<pt::id_t> s;
    for (auto [k,_] : this->data_ ) {
      s.insert(k);
    }
    return s;
  }

  const VtxRefVisits &get_ref_steps(pt::id_t ref_id) const {
    if (!pv_cmp::contains(this->data_, ref_id)) {
      // return an empty vector if the ref_id is not found
      static const VtxRefVisits empty_visits;
      return empty_visits;
    }
    return this->data_.at(ref_id);
  }

  // ---------
  // setter(s)
  // ---------

  void add(pt::id_t ref_id, pvt::AS s) {
    auto &visits = data_[ref_id]; // grab the reference once

    // append based on locus
    // find the first element whose step_idx >= s.step_idx()
    auto it = std::lower_bound(
        visits.begin(), visits.end(), s.get_step_idx(),
        [](auto const &a, auto step) { return a.get_step_idx() < step; });

    // insert before that position (or at end if none is >=)
    visits.insert(it, std::move(s));
  }
};

/**
  * A class to store the metadata on ref visits for a walk.
  *
  * It contains a set of ref ids that are present in the walk and a vector of
  * maps where each map is a step in the walk, where the index in the vector is
  * the step idx.
  *
  * Each map is a VtxRefMeta and thus contains a key of ref id and a value is
  * VtxRefVisits
 */
class WalkRefMeta {
  // the set of ref ids that are present in the walk
  std::set<id_t> ref_ids_;

  // a vector of maps where each map is a step in the walk and the index is the
  // step idx.The map contains a key of ref id and a value of a vector of
  // allele steps(AS) that the ref takes in the vertex at that step
  std::vector<VtxRefMeta> walk_ref_maps_;

public:
  // --------------
  // constructor(s)
  // --------------

  WalkRefMeta() {}

  // ---------
  // getter(s)
  // ---------

  const std::set<id_t> &get_ref_ids() const { return this->ref_ids_; }

  const VtxRefVisits &get_step_ref_data(pt::id_t ref_id, pt::idx_t step_idx) const {
    const VtxRefMeta &vrm = this->walk_ref_maps_.at(step_idx);
    return vrm.get_ref_steps(ref_id);
  }

  // ---------
  // setter(s)
  // ---------

  void add_vtx_ref_meta(VtxRefMeta &&vrm) {
     // extend ref_ids_ with the refs in the vrm
    for (const auto &ref_id : vrm.get_refs()) {
      this->ref_ids_.insert(ref_id);
    }
    this->walk_ref_maps_.push_back(vrm);
  }

  void add_ref_id(pt::id_t ref_id) { this->ref_ids_.insert(ref_id); }
};


/**
 * for a given AS in a walk, find the associated vertex in the graph and sort
 * associate the ref id and the steps, sort the steps in ascending order.
 */
VtxRefMeta get_vtx_itn(const bd::VG &g, const pvt::step_t &s) {
  pt::id_t v_id = s.v_id;
  const bd::Vertex &v = g.get_vertex_by_id(v_id);
  std::vector<bd::RefInfo> v_ref_data = v.get_refs();

  VtxRefMeta vrm;
  for (const bd::RefInfo &ref : v_ref_data) {
    pt::id_t ref_id = ref.get_ref_id();
    vrm.add(ref_id, pvt::AS::given_ref_info(v_id, ref));
  }

  return vrm;
}


/**
 * returns a pair of:
 * 1. the set of ref ids that are present in the walk
 * 2. a vector of maps where each map is a step in the walk and the index is the
 * step idx. The map contains a key of ref id and a value of a vector of
 * allele steps (AS) that the ref takes in the vertex at that step.
*/
WalkRefMeta comp_walk_ref_meta(const bd::VG &g, const pvt::walk_t &w) {
  WalkRefMeta wrm;

  for (pt::idx_t step_idx{}; step_idx < w.size(); ++step_idx) {
    const pvt::step_t &s = w[step_idx];
    VtxRefMeta ref_map = get_vtx_itn(g, s);
    wrm.add_vtx_ref_meta(std::move(ref_map));
  }

  return wrm;
}


pt::idx_t get_vtx_len(const bd::VG &g, pvt::AS s) {
  pt::id_t v_id = s.get_v_id();
  const bd::Vertex &prev_v = g.get_vertex_by_id(v_id);
  return prev_v.get_label().length();
}


/**
 * @brief compute min_locus and loop_no
 */
std::pair<pt::idx_t, pt::idx_t> comp_ref_visit_bounds(pt::id_t ref_id,
                                                      const WalkRefMeta &wrm,
                                                      const pvt::walk_t &w) {
  std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  // TODO: paralleise

  // minimum locus for the ref in the walk
  pt::idx_t min_locus{pc::MAX_IDX};
  for (pt::idx_t step_idx{}; step_idx < w.size(); ++step_idx) {

    const VtxRefVisits &vtx_ref_visits = wrm.get_step_ref_data(ref_id, step_idx);

    if (vtx_ref_visits.empty()) { // no visits for this ref in this step
      continue;
    }

    // because vtx_ref_visits is sorted by step_idx we can use the first element
    if (vtx_ref_visits.front().get_step_idx() < min_locus) {
      // ideally we can stop at the first find because it is the minimum because
      // the ref visits are sorted by locus
      min_locus = vtx_ref_visits.front().get_step_idx();
      break;
    }
  }

  // the number of times the ref is seen in the walk
  // this is also the step with the max ref visits
  pt::idx_t loop_no{0};

  for (pt::idx_t step_idx{}; step_idx < w.size(); ++step_idx) {

    const VtxRefVisits &vtx_ref_visits = wrm.get_step_ref_data(ref_id, step_idx);

    if (vtx_ref_visits.empty()) {
      // no visits for this ref in this step
      continue;
    }

    if (vtx_ref_visits.size() > loop_no) {
      loop_no = vtx_ref_visits.size();
    }
  }

  return {min_locus, loop_no};
}


/**
  * find the start locus for /subsequent/ loops that is greater than first loop
  * find curr locus for next loop
  * the next loop must start at the first step in the walk
  */
pt::idx_t find_loop_start_locus(const WalkRefMeta &wrm, pt::id_t ref_id, pt::idx_t curr_locus) {
  std::string fn_name{pv_cmp::format("[{}::{}]", MODULE, __func__)};

  const VtxRefVisits &vtx_ref_visits = wrm.get_step_ref_data(ref_id, 0);
  // get lowest number that is greater than curr_locus
  auto it =
    std::lower_bound(vtx_ref_visits.begin(), vtx_ref_visits.end(),
                     curr_locus, [](const pvt::AS &a, pt::idx_t locus) {
                       return a.get_step_idx() < locus;
                     });

  if (it == vtx_ref_visits.end()) {
    // throw an exception because we broke the assumption that the ref
    //  must start at the first step in other loops
    std::string msg =
      pv_cmp::format("{}: could not find next locus for ref {} in rov: ", fn_name, ref_id);
    throw std::runtime_error(msg);
  }

  return it->get_step_idx();
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

  WalkRefMeta wrm = comp_walk_ref_meta(g, w);

  for (pt::id_t ref_id : wrm.get_ref_ids()) { // for each ref in the walk

    auto [curr_locus, loop_count] = comp_ref_visit_bounds(ref_id, wrm, w);

    for (pt::idx_t loop_no{}; loop_no < loop_count; loop_no++) {
      pvt::AW allele_walk{w_idx};

      bool is_ref_cont{false}; // is the ref continuous in the walk?

      for (pt::idx_t step_idx{}; step_idx < w.size(); ++step_idx) {

        is_ref_cont = false; // reset for each step

        const VtxRefVisits &vtx_ref_visits = wrm.get_step_ref_data(ref_id, step_idx);

        for (const auto &s : vtx_ref_visits) {
          if (curr_locus == s.get_step_idx()) {
            allele_walk.append_step(s);
            curr_locus += get_vtx_len(g, s);
            is_ref_cont = true; // we have a step for this ref in the walk
          }
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
        curr_locus = find_loop_start_locus(wrm, ref_id, curr_locus);
      }
    }
  }

  return;
}

} // namespace variants

namespace graph {

// Maximum number of steps to take from flubble start to end
const pt::idx_t MAX_FLUBBLE_STEPS{20};

// direction for traversing a vertex in a bidirected graph
enum class dir_e {
  in,
  out
};

const dir_e IN = dir_e::in;
const dir_e OUT = dir_e::out;

typedef id_or_t idx_or_t; // specifically for idx instead of id

/**
  *@brief get the edges of a vertex end
  *
  *@param v the vertex
  *@param ve the vertex end
  *@return a set of edge indices
  */
inline const std::set<pt::idx_t> edges_at_end(const bd::Vertex &v, pgt::v_end_e ve) noexcept {
  return ve == pgt::v_end_e::l ? v.get_edges_l() : v.get_edges_r();
}

/**
  *@brief get vertex end from the traversal orientation and direction
  *
  *@param o the orientation
  *@param e the direction
  *@return the end of the vertex
  */
inline v_end_e get_v_end(or_e o, dir_e e) {
  switch (e) {
  case dir_e::in:
    return (o == or_e::forward) ? v_end_e::l : v_end_e::r;
  case dir_e::out:
    return (o == or_e::forward) ? v_end_e::r : v_end_e::l;
  default:
    std::string msg = pv_cmp::format("[{}::{}] Invalid v_end: {}",
                                  MODULE, __func__, static_cast<int>(e));
    throw std::invalid_argument(msg);
  }
};


/**
 *@brief get the orientation of a vertex end based on the side and direction of traversal
 *
 *@param side the side of the vertex end
 *@param d the direction of traversal
 *@return the orientation of the vertex end
 */
inline or_e get_or(pgt::v_end_e side, dir_e d) {
  switch (d) {
  case IN:
    return (side == pgt::v_end_e::l ? pgt::or_e::forward : pgt::or_e::reverse);
  case OUT:
    return (side == pgt::v_end_e::l ? pgt::or_e::reverse : pgt::or_e::forward);
  default:
    std::string msg = pv_cmp::format("[{}::{}] Invalid orientation: {}",
                  MODULE, __func__, static_cast<int>(d));
    throw std::invalid_argument(msg);
  }
};


/**
 * @brief the stack is a unique path from s to t
 */
pvt::walk_t walk_from_stack(const bd::VG &g, const std::deque<idx_or_t> &dq) {
  const std::string fn_name = pv_cmp::format("[povu::bidirected::{}]", __func__);

  //pvt::AW w;
  pvt::walk_t w;
  for (auto it = dq.begin(); it != dq.end(); ++it) {
    auto [v_idx, o] = *it;
    pgt::id_or_t s {g.v_idx_to_id(v_idx), o};
    // pvt::AS s{g.v_idx_to_id(v_idx), o};
    w.push_back(s);
  }

  return w;
}

// compute walks for flubble like
// Modified Johnson B algo to find paths between a pair of nodes in a graph
// Enumerate all simple s→t paths in a directed graph (adjacency lists).
// Runs in O((V+E)·(P+1)) time where P = number of paths found.
// when a vertex has been explored we unblock its neighbours and the current vertex
void comp_walks_fl_like(const bd::VG &g, pvt::RoV &rov) {

  std::vector<pvt::walk_t> &walks = rov.get_walks_mut();
  pvst::traversal_params_t tp = rov.get_pvst_vtx()->get_traversal_params();
  std::string label = rov.get_pvst_vtx()->as_str();

  auto [start_id, start_o] = tp.start;
  auto [stop_id, stop_o] = tp.end;

  pt::idx_t start_idx = g.v_id_to_idx(start_id);
  pt::idx_t stop_idx = g.v_id_to_idx(stop_id);

  idx_or_t s = { start_idx, start_o };
  idx_or_t t = { stop_idx, stop_o };

  std::deque<idx_or_t> dq;

  // the key is a vertex and the value is a set of vertices that have been seen
  // from that vertex
  std::map<id_or_t, std::set<idx_or_t>> seen_;

  idx_or_t curr = s;
  dq.push_back(curr);

  while (!dq.empty()) {
    // get the incoming vertices based on orientation
    curr = dq.back();

    if (curr == t) {
      walks.push_back(walk_from_stack(g, dq));
      dq.pop_back(); // we are done with this path up to this point
      continue; // we are done with that path
    }

    if (dq.size() > MAX_FLUBBLE_STEPS) {
      WARN("max steps reached for {}\n", rov.get_pvst_vtx()->as_str());
      return;
    }

    auto [v_idx, o] = curr;

    // if we have explored all neighbours of the current vertex
    bool is_explored {true};

    pgt::v_end_e ve = get_v_end(o, OUT);
    const bd::Vertex &v = g.get_vertex_by_idx(v_idx);
    const std::set<pt::idx_t> &nbr_edges = edges_at_end(v, ve);

    for (pt::idx_t e_idx : nbr_edges) {
      const bd::Edge &e = g.get_edge(e_idx);
      auto [side, alt_idx] = e.get_other_vtx(v_idx, ve);
      idx_or_t nbr {alt_idx, get_or(side, IN)};

      if(pv_cmp::contains(seen_[curr], nbr)) {
        continue;
      }

      is_explored = false;
      dq.push_back(nbr);
      seen_[curr].insert(nbr);
      break;
    }

    if (is_explored) {
      // unblock all in the stack up to and including the current vertex
      for (auto it = dq.rbegin(); it != dq.rend(); ++it) {
        seen_[*it].clear();
        if (*it == curr) {
          break;
        }
      }

      dq.pop_back(); // we are done with this path up to this point
    }
  }

  return;
}

void find_walks(const bd::VG &g, pvt::RoV &rov) {
  graph::comp_walks_fl_like(g, rov);
}
} // namespace graph
} // namespace povu::graph_utils
