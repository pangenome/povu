#include "./flubbles.hpp"
#include <cassert>
#include <utility>
#include <vector>

namespace povu::flubbles {

// orientation, id, class
struct oic_t {
  pgt::or_e orientation;
  pt::id_t id;
  pt::id_t st_idx;
  pt::idx_t cls;
};


struct eq_class_stack_t {
  std::vector<oic_t> s; // stack of equivalence classes
  std::vector<pt::idx_t> next_seen; // next seen index for each equivalence class

  // constructor
  eq_class_stack_t(pt::idx_t exp_size) {
    s.reserve(exp_size); // reserve some space for the stack
    next_seen.reserve(exp_size); // reserve some space for the next seen vector
  }
};

inline pvst::Vertex gen_fl(pt::id_t start_id, pgt::or_e start_or, pt::id_t end_id,
                               pgt::or_e end_or, pt::idx_t ai, pt::idx_t zi) {
  std::string fn_name = std::format("[povu::algorithms::flubble_tree::{}]", __func__);

  auto [a_g_id, a_or, z_g_id, z_or] =
      (start_or == pgt::or_e::reverse && end_or == pgt::or_e::reverse)
          ? std::make_tuple(end_id, pgt::or_e::forward, start_id, pgt::or_e::forward)
          : std::make_tuple(start_id, start_or, end_id, end_or);

  pgt::id_or_t a{a_g_id, a_or};
  pgt::id_or_t z{z_g_id, z_or};

  return pvst::Vertex::make_flubble(a, z, ai, zi);
}

/**
  * @brief Compute the ai and zi indices for a flubble
  *
  * @param st The spanning tree
  * @param a_idx The index of the first vertex (a)
  * @param z_idx The index of the second vertex (z)
  * @return A pair of indices (ai, zi) where ai is the index of the ancestor
  *         and zi is the index of the descendant in the spanning tree
 */
std::pair<pt::idx_t, pt::idx_t> compute_ai_zi(const pst::Tree &st,
                                              pt::idx_t a_idx,
                                              pt::idx_t z_idx) {

  std::string fn_name = std::format("[povu::algorithms::flubble_tree::{}]", __func__);

  std::vector<pt::idx_t> vtxs {a_idx, z_idx};

  auto get_vtx_pair = [&](pt::idx_t v_idx) -> void {
    if (!st.is_root(v_idx)) {
      const pst::Edge &e = st.get_parent_edge(v_idx);
      if (e.get_color() == pgt::color_e::black) {
        pt::idx_t v_idx_ = st.get_parent_v_idx(v_idx);
        vtxs.push_back(v_idx_);
        return;
      }
      else {
        for (auto c_e_idx : st.get_child_edge_idxs(v_idx)) {
          const pst::Edge &ce = st.get_tree_edge(c_e_idx);
          if (ce.get_color() == pgt::color_e::black) {
            pt::idx_t v_idx_ = ce.get_child_v_idx();
            vtxs.push_back(v_idx_);
            return;
          }
        }
      }
    }
  };

  get_vtx_pair(a_idx);
  get_vtx_pair(z_idx);
  std::sort(vtxs.begin(), vtxs.end(), [&](pt::idx_t a, pt::idx_t b) { return a < b; });

#ifdef DEBUG // vtx should contain exactly 4 elements
  assert(vtxs.size() == 4);
#endif

  return std::make_pair(vtxs[1], vtxs[2]);
}

/**
  * @brief
 */
void add_flubbles(const pst::Tree &st, const eq_class_stack_t &ecs,
                  pvtr::Tree<pvst::Vertex> &vst) {
  std::string fn_name = std::format("[povu::algorithms::flubble_tree::{}]", __func__);

  const auto &[stack_, next_seen] = ecs;

  struct ci {
    pt::idx_t cl; // class
    pt::idx_t idx; // index in stack_
  };

  std::stack<ci> s;
  std::unordered_set<pt::idx_t> in_s; // classes in s

  pt::idx_t prt_v { vst.root_idx() }; // parent vertex

  //pt::idx_t ctr = vst.vtx_count();

  for (pt::idx_t i {}; i < stack_.size(); ++i) {

    auto [or_curr, id_curr, st_idx_curr, cl_curr] = stack_[i];

    if (id_curr == pc::INVALID_IDX) {
      continue;
    }

    // find the parent vertex, applies for non-siblings
    if (in_s.contains(cl_curr)) {

      // pop until (and including) the one whose cl equals cl_curr
      while(!s.empty()) {
        auto [cl, _] = s.top();
        s.pop();
        in_s.erase(cl);

        if (cl == cl_curr){
          break;
        }
      }

      if (prt_v != vst.root_idx()) {
        prt_v = vst.get_parent_idx(prt_v);
      }
    }

    if ((i + 1) < next_seen[i]) {
      auto [or_nxt, id_nxt, st_idx_nxt, _] = stack_[next_seen[i]];

      if (id_nxt == pc::INVALID_IDX) {
        continue;
      }

      auto [ai, zi] = compute_ai_zi(st, st_idx_curr, st_idx_nxt);

      pvst::Vertex vtx = gen_fl(id_curr, or_curr, id_nxt, or_nxt, ai, zi);
      pt::idx_t pvst_v_idx = vst.add_vertex(vtx);
      vst.add_edge(prt_v, pvst_v_idx);
      prt_v = pvst_v_idx;
    }

    s.push( {cl_curr, i} );
    in_s.insert(cl_curr);
  }
}

/**
 * @brief Compute the next seen index for each equivalence class
 *
 * @param stack_ The equivalence class stack
 * @return The next seen index for each equivalence class
 */
void compute_eq_class_metadata(eq_class_stack_t &ecs) {
  std::string fn_name = std::format("[povu::algorithms::flubble_tree::{}]", __func__);

  const std::vector<oic_t> &stack_ = ecs.s;
  std::vector<pt::idx_t> &next_seen = ecs.next_seen;

  for (pt::idx_t i {}; i < stack_.size(); ++i) {
    next_seen.push_back(pc::INVALID_CLS);
  }

  // an eq class and the index of the next time it is encountered in stack_
  std::unordered_map<pt::idx_t, pt::idx_t> next_seen_map;
  next_seen_map.reserve(stack_.size());

  for (std::size_t i {stack_.size()};  i-- > 0; ) {
    auto [or_curr, id_curr, _, cl_curr] = stack_[i];
    pt::idx_t next_idx = next_seen_map.contains(cl_curr) ? next_seen_map[cl_curr] : i;

    next_seen[i] = next_idx;
    next_seen_map[cl_curr] = i;
  }
}


void compute_eq_class_stack(const pst::Tree &st, std::vector<oic_t> &stack) {
  std::string fn_name = std::format("[povu::algorithms::flubble_tree::{}]", __func__);

  ptu::BranchDesc desc = ptu::br_desc(st);

  auto is_branching = [&](pt::idx_t v_idx) -> bool {
    return (st.get_children(v_idx).size() > 1);
  };

  pt::idx_t root_idx = st.get_root_idx();

  // TODO: use desc to get these
  // edge idx => stackette
  typedef std::map<pt::idx_t, std::list<oic_t>> edge_stack_t ;
  // v idx to mini stack
  std::map<pt::idx_t, edge_stack_t> cache; // edge_map

  std::list<oic_t> mini_stack; // TODO: rename

  auto merge_stacks = [&](std::list<oic_t> &from, std::list<oic_t> &into) {
    // put from infront of or on top of into
    into.splice(into.begin(), from);
  };

  for (pt::idx_t v_idx{st.vtx_count()}; v_idx-- > 0; ) {

    if (v_idx == root_idx || is_branching(v_idx)) {

      // merge from cache into stack_map
      edge_stack_t &stackettes = cache[v_idx];
      const std::vector<pt::idx_t> &c_edges = desc[v_idx].sorted_br;

      for (auto e_idx : c_edges) {
        std::list<oic_t> &stackette = stackettes[e_idx];
        merge_stacks(stackette, mini_stack);
      }
    }

    if (v_idx == root_idx) {
      break; // done
    }

    const pst::Edge &e = st.get_parent_edge(v_idx);

    if (e.get_color() == color::black) {
      or_e o = st.get_vertex(v_idx).type() == pgt::v_type_e::r ? pgt::or_e::forward : pgt::or_e::reverse;
      mini_stack.push_front({o, st.get_vertex(v_idx).g_v_id(), v_idx, e.get_class()});
    }

    // if the parent is braching
    if (is_branching(st.get_parent_v_idx(v_idx))) {
      pt::idx_t e_idx = st.get_vertex(v_idx).get_parent_e_idx();
      pt::idx_t p_v_idx = st.get_parent_v_idx(v_idx);

      edge_stack_t &es = (cache.contains(p_v_idx)) ? cache[p_v_idx] : cache[p_v_idx] = {};

      es[e_idx] = mini_stack;
      cache[p_v_idx] = es;
      mini_stack.clear();
    }
  }

  // convert the mini_stack to stack
  for (auto it = mini_stack.begin(); it != mini_stack.end(); ++it) {
    stack.emplace_back(*it);
  }
}


pvtr::Tree<pvst::Vertex> find_flubbles(const pst::Tree &st) {
  std::string fn_name = std::format("[povu::algorithms::{}]", __func__);

  eq_class_stack_t ecs{ st.tree_edge_count() };
  {
    compute_eq_class_stack(st, ecs.s);
    compute_eq_class_metadata(ecs);
  }

  // create the pvst and add the root vertex
  pvtr::Tree<pvst::Vertex> vst;
  {
    pvst::Vertex root_v = pvst::Vertex::make_dummy();
    pt::idx_t root_v_idx = vst.add_vertex(root_v);
    vst.set_root_idx(root_v_idx);
  }
  add_flubbles(st, ecs, vst);

  return vst;
}
} // namespace povu::flubbles
