#include "./flubbles.hpp"
#include <cassert>
#include <utility>
#include <vector>

namespace povu::flubbles {

pvst::Flubble gen_fl(pt::id_t start_id, pgt::or_e start_or, pt::id_t end_id,
                               pgt::or_e end_or, pt::idx_t ai, pt::idx_t zi) {
  std::string fn_name = std::format("[povu::algorithms::flubble_tree::{}]", __func__);

  auto [a_g_id, a_or, z_g_id, z_or] =
      (start_or == pgt::or_e::reverse && end_or == pgt::or_e::reverse)
          ? std::make_tuple(end_id, pgt::or_e::forward, start_id, pgt::or_e::forward)
          : std::make_tuple(start_id, start_or, end_id, end_or);

  pgt::id_or_t a{a_g_id, a_or};
  pgt::id_or_t z{z_g_id, z_or};

  return pvst::Flubble(a, z, ai, zi);
  //return pvst::Vertex::make_flubble(a, z, ai, zi);
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
std::pair<pt::idx_t, pt::idx_t> compute_ai_zi(const pst::Tree &st, pt::idx_t a_e_idx, pt::idx_t z_e_idx) {

  std::string fn_name = std::format("[povu::algorithms::flubble_tree::{}]", __func__);

  std::vector<pt::idx_t> vtxs {};
  auto get_vtx_pair = [&](pt::idx_t e_idx) -> void {
    const pst::Edge &e = st.get_tree_edge(e_idx);
    vtxs.push_back(e.get_child_v_idx());
    vtxs.push_back(e.get_parent_v_idx());
  };

  get_vtx_pair(a_e_idx);
  get_vtx_pair(z_e_idx);

  // sort the vertices by their index in the spanning tree
  std::sort(vtxs.begin(), vtxs.end(), [&](pt::idx_t a, pt::idx_t b) { return a < b; });

#ifdef DEBUG // vtx should contain exactly 4 elements
  assert(vtxs.size() == 4);
#endif

  if (vtxs[1] == 1096 && vtxs[2] == 1099) {
    std::cerr << fn_name << "..." << "\n";

    for (auto v : vtxs) {
      std::cerr << v << "(" << st.get_vertex(v).g_v_id() << "). Children: ";
      // print the children
      for (auto c : st.get_children(v)) {
        std::cerr << c << "(" << st.get_vertex(c).g_v_id() << "), ";
      }
      std::cerr << "\n";
    }

    

    std::cerr << "\n";
  }


  return std::make_pair(vtxs[1], vtxs[2]);
}

/**
  * @brief
 */
void add_flubbles(const pst::Tree &st, const eq_class_stack_t &ecs, pvtr::Tree &vst) {
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

      pvst::Flubble vtx = gen_fl(id_curr, or_curr, id_nxt, or_nxt, ai, zi);
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
    pt::idx_t pe = st.get_vertex(v_idx).get_parent_e_idx();

    if (e.get_color() == color::black) {
      or_e o = st.get_vertex(v_idx).type() == pgt::v_type_e::r ? pgt::or_e::forward : pgt::or_e::reverse;
      mini_stack.push_front({o, st.get_vertex(v_idx).g_v_id(), pe, e.get_class()});
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


void handle_vertex(pst::Tree &t,
                   std::size_t v,
                   std::vector<boundary> &hairpins,
                   boundary &curr_bry,
                   bool &in_hairpin,
                   std::set<std::size_t> &articulated_vertices) {
  std::string fn_name = std::format("[povu::algorithms::{}]", __func__);

  /*
   * compute v.hi
   * ------------
   */

  pt::idx_t hi_0 {pc::INVALID_IDX};
  std::set<std::size_t> obe = t.get_obe(v);
  for (auto be : obe) {
    hi_0 = std::min(hi_0, t.get_vertex(be).dfs_num());
  }

  // given a node v find its child with the lowest hi value
  // (closest to root)
  // its hi value is hi_1 and the dfs num of that vertex is hi_child
  // children are empty for dummy stop node

  pt::idx_t hi_1 { pc::INVALID_IDX };
  std::set<std::size_t> children = t.get_children(v);

  bool is_leaf = children.empty();
  // insert current boundary into the boundary list
  // if we are in a hairpin and
  // if we are in a leaf or got to the root
  if (in_hairpin && ((is_leaf && !t.is_root(v)) || t.is_root(v))) {
    hairpins.push_back(std::move(curr_bry));
    curr_bry= NULL_BOUNDARY;
    in_hairpin = false;
  }


  // a vector of pairs of hi values and children
  std::vector<std::pair<std::size_t, std::size_t>> hi_and_child{};
  hi_and_child.reserve(children.size());

  // for each child input the hi value and the child
  for (std::size_t child: children) {
    hi_and_child.push_back({t.get_vertex(child).hi(), child});
  }
  std::sort(hi_and_child.begin(), hi_and_child.end());
  if (!children.empty()) { hi_1 = hi_and_child[0].first; }

  t.get_vertex_mut(v).set_hi(std::min(hi_0, hi_1));

  // the child vertex whose hi value is hi_1
  std::size_t hi_child { pc::INVALID_IDX };
  for (std::size_t child: children) {
    if (t.get_vertex(child).hi() == hi_1) {
      hi_child = child;
      break;
    }
  }

  // if hi_and_child has at least 2 elements
  // assign hi_2 from the second element of hi_and_child
  // works because hi_and_child is sorted
  std::size_t hi_2 { pc::INVALID_IDX };
  for (std::size_t child: children) {
    if (child != hi_child && t.get_vertex(child).hi() < v && articulated_vertices.find(child) == articulated_vertices.end()) {
      hi_2 = t.get_vertex(child).hi();
      break;
    }
  }

  //std::cout << "\tcompute bracket list";
  /*
   * compute bracket list
   * --------------------
   */

  // the bracketlist was created in tree constructor
  //if (children.size() > 1) {
  //std::cerr << "many ~> " << children.size() << std::endl;
  // }
  for (auto c: children) {
    t.concat_bracket_lists(v, c);
  }



  // for each capping backedge TODO: add

  // pop incoming backedges
  // remove backedges we have reached the end of
  std::set<std::size_t> ibe_idxs = t.get_ibe_idxs(v);
  for (std::size_t b: ibe_idxs) {
    t.del_bracket(v, b);

    // TODO: set backedge class ?? was id not enough?
    // do this in the del_bracket_method?
    pst::BackEdge& be= t.get_backedge(b);
    if (be.type() != pst::be_type_e::capping_back_edge && !be.is_class_defined()) {
      be.set_class(t.new_class());
    }
  }


  // push outgoing backedges
  std::set<std::size_t> obe_i = t.get_obe_idxs(v);
  for (std::size_t be_idx : obe_i) {
    t.push(v, be_idx);
  }


  if (hi_2 < hi_0) {
    // add a capping backedge
    std::size_t dest_v =  hi_2;
    std::size_t be_idx = t.add_be(v, dest_v, pst::be_type_e::capping_back_edge, color_e::gray);
    t.push(v, be_idx);
  }


  if (t.get_bracket_list(v).empty()) {
    std::size_t dest_v = t.get_root_idx();
    if (t.get_vertex(v).type() != v_type_e::dummy) {
      //std::cerr << "add art be " << t.get_vertex(v).name() << " " << dest_v << std::endl;
      if (curr_bry.b1 != pc::INVALID_IDX) { std::cerr << fn_name << "WARN: curr boundary already set\n"; }
      curr_bry.b1 = t.get_vertex(v).g_v_id();
      //std::cerr << "Found hairpin boundary start " << t.get_vertex(v).g_v_id() << std::endl;
    }

    // add a simplifying back edge
    std::size_t be_idx = t.add_be(v, dest_v, pst::be_type_e::simplifying_back_edge, color_e::gray);
    t.push(v, be_idx);
    t.get_vertex_mut(v).set_hi(t.get_root_idx());

    in_hairpin = true;
  }
  else if (in_hairpin) {

    pst::Bracket& b = t.top(v);

    std::size_t b_id = b.back_edge_id();
    pst::BackEdge &be = t.get_backedge_ref_given_id(b_id);

    // extend the end boudary of the current hairpin
    if (be.type() == pst::be_type_e::simplifying_back_edge) {
      //boundary = v;
      curr_bry.b2 = t.get_vertex(v).g_v_id();
    }

  }

  /*
   * determine equivalance class for edge v.parent() to v
   * ---------------------------------------------------
   */

  // if v is not the root of the spanning tree
  if (!t.is_root(v)) {

    /*default behavior*/

    pst::Bracket& b = t.top(v);

    if (t.list_size(v) !=  b.recent_size()) {
      b.set_recent_size(t.list_size(v));
      b.set_recent_class(t.new_class());
    }

    // when retreating out of a node the tree edge is labelled with
    // the class of the topmost bracket in the bracket stack
    pst::Edge& e = t.get_incoming_edge(v);
    e.set_class(b.recent_class());

    /*check for e, b equivalance*/
    if (b.recent_size() == 1) {
      std::size_t b_id = b.back_edge_id();
      pst::BackEdge& be = t.get_backedge_ref_given_id(b_id);
      be.set_class(e.get_class());
    }
  }
}

void simple_cycle_equiv(pst::Tree &t, const core::config &app_config) {

  std::string fn_name = std::format("[povu::algorithms::{}]", __func__);

  std::set<std::size_t> articulated_vertices;

  std::vector<boundary> boundaries;
  boundaries.reserve(EXPECTED_HAIRPIN_COUNT);
  bool in_hairpin { false };
  boundary curr_bry { NULL_BOUNDARY }; // current_boundary


  for (pt::idx_t v {t.vtx_count() - 1}; v < pc::MAX_IDX; --v) {
    handle_vertex(t, v, boundaries, curr_bry, in_hairpin, articulated_vertices);
  }

  // print boundaries
  if (app_config.inc_hairpins()) {
    for (auto b : boundaries) {
      std::cerr << "Boundary: " << b.b1 << " " << b.b2 << std::endl;
    }
  }

}

pvtr::Tree find_flubbles(pst::Tree &st, const core::config &app_config) {
  std::string fn_name = std::format("[povu::algorithms::{}]", __func__);

  simple_cycle_equiv(st, app_config);

  eq_class_stack_t ecs{st.tree_edge_count()};
  {
    compute_eq_class_stack(st, ecs.s);
    compute_eq_class_metadata(ecs);
  }

  // create the pvst and add the root vertex
  pvtr::Tree vst;
  {
    pvst::Dummy root_v;
    pt::idx_t root_v_idx = vst.add_vertex(root_v);
    vst.set_root_idx(root_v_idx);
  }
  add_flubbles(st, ecs, vst);

  return vst;
}

} // namespace povu::flubbles
