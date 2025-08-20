#include "./graph.hpp"
#include <utility>

namespace povu::genomics::graph {

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
inline pgt::v_end_e get_v_end(pgt::or_e o, dir_e e) {
  switch (e) {
  case dir_e::in:
    return (o == pgt::or_e::forward) ? pgt::v_end_e::l : pgt::v_end_e::r;
  case dir_e::out:
    return (o == pgt::or_e::forward) ? pgt::v_end_e::r : pgt::v_end_e::l;
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
inline pgt::or_e get_or(pgt::v_end_e side, dir_e d) {
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
pvt::walk_t walk_from_stack(const bd::VG &g, const std::deque<idx_or_t> &dq,
                            pvst::route_e route) {
  pvt::walk_t w;
  auto append_w = [&](idx_or_t it) {
    auto [v_idx, o] = it;
    pgt::id_or_t s{g.v_idx_to_id(v_idx), o};
    w.push_back(s);
  };

  switch (route) {
  case pvst::route_e::s2e:
    for (auto &&it : dq) append_w(it);
    break;
  case pvst::route_e::e2s:
    for (auto it = dq.crbegin(); it != dq.crend(); ++it) append_w(*it);
    break;
  default:
    ERR("Unsupported route type: {}\n", static_cast<int>(route));
    std::exit(1);
  }

  return w;
}

// compute walks for flubble like
// Modified Johnson B algo to find paths between a pair of nodes in a graph
// Enumerate all simple s→t paths in a directed graph (adjacency lists).
// Runs in O((V+E)·(P+1)) time where P = number of paths found.
// when a vertex has been explored we unblock its neighbours and the current vertex
void comp_walks_s2e(const bd::VG &g, idx_or_t s, idx_or_t t,
                    std::vector<pvt::walk_t> &walks,
                    const std::string_view &rov_label) {

  std::deque<idx_or_t> dq;

  // the key is a vertex and the value is a set of vertices that have been seen
  // from that vertex
  std::map<pgt::id_or_t, std::set<idx_or_t>> seen_;

  idx_or_t curr = s;
  dq.push_back(curr);

  while (!dq.empty()) {
    // get the incoming vertices based on orientation
    curr = dq.back();

    if (curr == t) {
      walks.push_back(walk_from_stack(g, dq, pvst::route_e::s2e));
      dq.pop_back(); // we are done with this path up to this point
      continue; // we are done with that path
    }

    if (dq.size() > MAX_FLUBBLE_STEPS) {
      WARN("max steps reached for {}\n", rov_label);
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

// compute walks for flubble like
// Modified Johnson B algo to find paths between a pair of nodes in a graph
// Enumerate all simple s→t paths in a directed graph (adjacency lists).
// Runs in O((V+E)·(P+1)) time where P = number of paths found.
// when a vertex has been explored we unblock its neighbours and the current vertex
void comp_walks_e2s(const bd::VG &g, idx_or_t s, idx_or_t t,
                    std::vector<pvt::walk_t> &walks, const std::string_view &rov_label) {

  std::deque<idx_or_t> dq;

  // the key is a vertex and the value is a set of vertices that have been seen
  // from that vertex
  std::map<pgt::id_or_t, std::set<idx_or_t>> seen_;

  idx_or_t curr = t;
  dq.push_back(curr);

  while (!dq.empty()) {
    // get the incoming vertices based on orientation
    curr = dq.back();

    if (curr == s) {
      walks.push_back(walk_from_stack(g, dq, pvst::route_e::e2s));
      dq.pop_back(); // we are done with this path up to this point
      continue; // we are done with that path
    }

    if (dq.size() > MAX_FLUBBLE_STEPS) {
      WARN("max steps reached for {}\n", rov_label);
      return;
    }

    auto [v_idx, o] = curr;

    // if we have explored all neighbours of the current vertex
    bool is_explored {true};

    pgt::v_end_e ve = get_v_end(o, IN);
    const bd::Vertex &v = g.get_vertex_by_idx(v_idx);
    const std::set<pt::idx_t> &nbr_edges = edges_at_end(v, ve);

    for (pt::idx_t e_idx : nbr_edges) {
      const bd::Edge &e = g.get_edge(e_idx);
      auto [side, alt_idx] = e.get_other_vtx(v_idx, ve);
      idx_or_t nbr {alt_idx, get_or(side, OUT)};

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
  // assume route params are set
  pvst::route_params_t rp = *rov.get_pvst_vtx()->get_route_params();
  auto [rp_l, rp_r, rp_route] = rp;
  auto [start_id, start_o] = rp_l;
  auto [stop_id, stop_o] = rp_r;

  pt::idx_t start_idx = g.v_id_to_idx(start_id);
  pt::idx_t stop_idx = g.v_id_to_idx(stop_id);

  idx_or_t s = {start_idx, start_o};
  idx_or_t t = {stop_idx, stop_o};

  switch (rp_route) {
  case pvst::route_e::s2e:
    return comp_walks_s2e(g, s, t, rov.get_walks_mut(), rov.get_pvst_vtx()->as_str());
  case pvst::route_e::e2s: // if you have a specific enum for this
    return comp_walks_e2s(g, s, t, rov.get_walks_mut(), rov.get_pvst_vtx()->as_str());
  default:
    ERR("Unsupported route type: {}\n", static_cast<int>(rp_route));
    std::exit(1);
  }
}

} // namespace povu::graph_utils
