#include "./graph_utils.hpp"
#include "types/graph.hpp"
#include "types/pvst.hpp"
#include <vector>

namespace povu::graph_utils {

/**
 *@brief a helper to generate a vector of walks
 *
 *@param preds a map of predecessors (value is all predecessors)
 *@param r a region of variation
 *@return a vector of walks
 */
std::vector<Walk> collapse(std::map<id_or_t, std::vector<id_or_t>> &preds, const pvst::Flubble *ft_v) {
  const std::string fn_name = std::format("[povu::bidirected::{}]", __func__);

  // auto [start_id, start_o] = r.get_entry();
  // auto [stop_id, stop_o] = r.get_exit();

  auto [start_id, start_o] = ft_v->get_a();
  auto [stop_id, stop_o] = ft_v->get_z();

  //std::vector<pgt::Walk> walks;

  // Each element in the stack holds a partial path (starting at the
  // destination).
  std::stack<pgt::Walk> stack;
  pgt::Walk w {stop_id, stop_o};
  stack.push(w);

  std::vector<pgt::Walk> allWalks;

  while (!stack.empty()) {
    pgt::Walk curr = stack.top();
    stack.pop();

    pgt::Step top_step = curr.get_steps().back();

    // std::cerr << "top " << top_step.get_o() << top_step.get_v_id() << "\n";

    if (top_step.get_v_id() == start_id && top_step.get_o() == start_o) {
      allWalks.push_back(curr);
    }

    auto preds_it = preds.find({top_step.get_v_id(), top_step.get_o()});
    if (preds_it == preds.end()) {
      //std::cerr << "not in map\n";
      continue;
    }
    else {
      //std::cerr << "in map\n";
    }

    for (auto pred : preds[{top_step.get_v_id(), top_step.get_o()}]) {
      auto [v_id, o] = pred;
      //std::cerr << "pred " << o << v_id << "\n";
      pgt::Step s {v_id, o};
      pgt::Walk w = curr;
      w.append_step(s);
      stack.push(w);
    }
  }

  for (auto &w: allWalks) {
    std::reverse(w.get_steps_mut().begin(), w.get_steps_mut().end());
  }

  // print the walks in all walks
  // for (auto &w: allWalks) {
  //   std::cerr << w.as_str() << "\n";
  // }

  return allWalks;
}

std::vector<pgt::Walk> populate_walks(const bd::VG &g,
                                      const pvst::VertexBase *pvst_vtx_ptr) {
  const std::string fn_name = std::format("[povu::bidirected::{}]", __func__);

  const pvst::Flubble *ft_v = static_cast<const pvst::Flubble *>(pvst_vtx_ptr);

  typedef id_or_t idx_or_t; // specifically for idx instead of id

   //  direction
  enum class dir_e {
    in,
    out
  };

  const dir_e IN = dir_e::in;
  const dir_e OUT = dir_e::out;

  auto [start_id, start_o] = ft_v->get_a();
  auto [stop_id, stop_o] = ft_v->get_z();

  pt::idx_t start_idx = g.v_id_to_idx(start_id);
  pt::idx_t stop_idx = g.v_id_to_idx(stop_id);

  idx_or_t s = {start_idx, start_o};
  idx_or_t t = {stop_idx, stop_o};

  std::queue<idx_or_t> q;

  // a map to keep track of the vertices whose incoming neighbours paths we have
  // extended so far key is the vertex and value is the set of vertices whose
  // paths we have extended
  std::map<idx_or_t, std::set<idx_or_t>> seen;

  // a set to keep track of the vertices we've seen
  std::set<idx_or_t> explored;

  //bool all_incoming_explored{true};

  //std::size_t counter {}; // a counter to keep track of the number of iterations
  // allows us to short circuit the traversal if counter > max_steps

  // a map of an idx and orientation to the incoming walk from the entry to the
  // idx and side that is incoming for the orientation
  // a predecessor map
  std::map<id_or_t, std::vector<id_or_t>> in_walks;

  auto get_v_end = [](or_e o, dir_e e) {
    return e == dir_e::in ? (o == or_e::forward ? v_end_e::l : v_end_e::r)
                          : (o == or_e::forward ? v_end_e::r : v_end_e::l);
  };


  auto get_alt_or = [](pgt::v_end_e side, dir_e d) -> pgt::or_e {
    return d == dir_e::in
      ? (side == pgt::v_end_e::r ? pgt::or_e::forward : pgt::or_e::reverse)
      : (side == pgt::v_end_e::l ? pgt::or_e::forward : pgt::or_e::reverse);
  };


  auto to_id_or = [&g](idx_or_t i) -> pgt::id_or_t {
    auto [v_idx, o] = i;
    return {g.v_idx_to_id(v_idx), o};
  };

  auto to_idx_or = [&g](pgt::id_or_t i) -> idx_or_t {
    auto [v_id, o] = i;
    return {g.v_id_to_idx(v_id), o};
  };

  idx_or_t curr = t;
  q.push(curr);

  //std::cerr << "handling: " << r.get_flb().as_str() << "\n";


  //std::cerr << "start " << s << "\n";

  while ( !q.empty() ) {
    // get the incoming vertices based on orientation
    idx_or_t curr = q.front();
    auto [v_idx, o] = curr;

    pgt::id_or_t curr_id_or = to_id_or(curr);
    q.pop();

    if (curr == s) {
      continue;
    }

    if (in_walks.size() > MAX_FLUBBLE_STEPS) {
      std::cerr << "max steps reached for " << pvst_vtx_ptr->as_str() << "\n";
      return std::vector<pgt::Walk>();
    }

    //std::cerr << "curr" << curr_id_or << std::endl;

    pgt::v_end_e ve = get_v_end(o, IN);
    const bd::Vertex &v = g.get_vertex_by_idx(v_idx);
    const std::set<pt::idx_t> &nbr_edges = ve == pgt::v_end_e::l ? v.get_edges_l() : v.get_edges_r();

    for (pt::idx_t e_idx : nbr_edges) {
      const bd::Edge &e = g.get_edge(e_idx);
      auto [side, alt_idx] = e.get_other_vtx(v_idx, ve);

      idx_or_t nbr {alt_idx, get_alt_or(side, IN)};

      if (seen[curr].contains(nbr)) {
        continue;
      }

      q.push(nbr);
      //std::cerr << "pushing " << to_id_or(nbr) << std::endl;

      in_walks[curr_id_or].push_back(to_id_or(nbr));

      seen[curr].insert(nbr);
    }
  }

  // std::cerr << "rov " << r.get_entry().as_str() << " in_walks size " << in_walks.size() << std::endl;

  //r.set_walks(collapse(in_walks, r));
  return collapse(in_walks, ft_v);
  //return;
}

std::vector<pgt::Walk> get_walks(const bd::VG &g, const pvst::VertexBase *pvst_vtx_ptr) {
  const std::string fn_name{std::format("[{}::{}]", MODULE, __func__)};

  std::vector<pgt::Walk> walks;

  if (pvst_vtx_ptr == nullptr) {
    throw std::invalid_argument(std::format("{}: pvst_vtx_ptr is null", fn_name));
  }

  pvst::traversal_params_t p = pvst_vtx_ptr->get_traversal_params();

  if (!p.traversable) {
    std::cerr << std::format("{}: vertex {} is not traversable\n", fn_name, pvst_vtx_ptr->as_str());
    return walks; // return empty vector if the vertex is not traversable
  }

  auto is_fl_like = [](pvst::vt_e typ) -> bool {
    return typ == pvst::vt_e::flubble || typ == pvst::vt_e::tiny || typ == pvst::vt_e::parallel;
  };

  if (is_fl_like(pvst_vtx_ptr->get_type())) {
    return populate_walks(g, pvst_vtx_ptr);
  }

  return walks;
}

} // namespace povu::graph_utils
