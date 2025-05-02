#include <algorithm>
#include <cstddef>
#include <iostream>

#include <unistd.h>
#include <utility>
#include <vector>
#include <format>
#include <chrono>

#include "./algorithms.hpp"


namespace povu::algorithms {

/**
 * Compute the equivalance class of a given vertex
 * Find the cycle equivalent edges
 * in reverse DFS
 */
void eulerian_cycle_equiv(pst::Tree &t) {
  std::string fn_name = std::format("[povu::algorithms::{}]", __func__);

  std::set<std::size_t> articulted_vertices;
  bool report_time {false};
  bool check_time{false};

  auto start = std::chrono::high_resolution_clock::now();

  bool in_hairpin {false};

  std::size_t boundary {pc::UNDEFINED_SIZE_T};

  for (std::size_t v { t.size() - 1 }; v < pc::UNDEFINED_SIZE_T; --v) {

    if (report_time && v % 10000 == 0) {
      check_time = true;
      std::cerr << "v: " << v << std::endl;
      start = std::chrono::high_resolution_clock::now();
    }


    // check if v is a multiple of 100000

    //t0 = pt::Time::now();

    /*
     * compute v.hi
     * ------------
     */

    pt::idx_t hi_0 { pc::UNDEFINED_IDX };
    std::set<std::size_t> obe = t.get_obe(v);
    for (auto be: obe) {
      hi_0 = std::min(hi_0, t.get_vertex(be).dfs_num());
    }

    // given a node v find its child with the lowest hi value
    // (closest to root)
    // its hi value is hi_1 and the dfs num of that vertex is hi_child
    // children are empty for dummy stop node
    pt::idx_t hi_1 { pc::UNDEFINED_IDX };

    std::set<std::size_t> children = t.get_children(v);


    if (in_hairpin && children.empty() && !t.is_root(v)) { // v is a leaf
      in_hairpin = false;
      std::cerr << "Found hairpin boundary end " << t.get_vertex(boundary).g_v_id() << std::endl;
    }
    else if (in_hairpin && t.is_root(v)) {
      in_hairpin = false;
      std::cerr << "Found hairpin boundary end " << t.get_vertex(boundary).g_v_id() << std::endl;
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
    std::size_t hi_child { pc::UNDEFINED_SIZE_T };
    for (std::size_t child: children) {
      if (t.get_vertex(child).hi() == hi_1) {
        hi_child = child;
        break;
      }
    }

    // if hi_and_child has at least 2 elements
    // assign hi_2 from the second element of hi_and_child
    // works because hi_and_child is sorted
    std::size_t hi_2 { pc::UNDEFINED_SIZE_T };
    for (std::size_t child: children) {
      if (child != hi_child && t.get_vertex(child).hi() < v && articulted_vertices.find(child) == articulted_vertices.end()) {
        hi_2 = t.get_vertex(child).hi();
        break;
      }
    }


    if (report_time && check_time) {
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::milli> duration1 = end - start;
      std::cout << "Block 1 Time: " << duration1.count() << " ms" << std::endl;

      start = std::chrono::high_resolution_clock::now();
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


    if (report_time && check_time) {
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::milli> duration1 = end - start;
      std::cout << "Block 2 concat Time: " << duration1.count() << " ms" << std::endl;

      start = std::chrono::high_resolution_clock::now();
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


    if (report_time && check_time) {
      auto end = std::chrono::high_resolution_clock::now();
      std::chrono::duration<double, std::milli> duration1 = end - start;
      std::cout << "Block 2 del Time: " << duration1.count() << " ms" << std::endl;

      start = std::chrono::high_resolution_clock::now();
    }


    // push outgoing backedges
    std::set<std::size_t> obe_i = t.get_obe_idxs(v);
    for (std::size_t be_idx : obe_i) {
      t.push(v, be_idx);
    }


    if (hi_2 < hi_0) {
      // add a capping backedge
      std::size_t dest_v =  hi_2;
      std::size_t be_idx =
        t.add_be(v, dest_v, pst::be_type_e::capping_back_edge, color_e::gray);
      t.push(v, be_idx);
    }



    if (t.get_bracket_list(v).empty()) {
      std::size_t dest_v = t.get_root_idx();
      if (t.get_vertex(v).type() != v_type_e::dummy) {
        //std::cerr << "add art be " << t.get_vertex(v).name() << " " << dest_v << std::endl;

        std::cerr << "Found hairpin boundary start " << t.get_vertex(v).g_v_id() << std::endl;
      }

      std::size_t be_idx = t.add_be(v, dest_v, pst::be_type_e::simplifying_back_edge, color_e::gray);
      t.push(v, be_idx);
      t.get_vertex_mut(v).set_hi(t.get_root_idx());

      in_hairpin = true;
    }
    else if (in_hairpin) {

      pst::Bracket& b = t.top(v);

      std::size_t b_id = b.back_edge_id();
      pst::BackEdge &be = t.get_backedge_ref_given_id(b_id);
      if (be.type() == pst::be_type_e::simplifying_back_edge) {
        boundary = v;
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

}

// a hairpin boundary
struct boundary {
  std::size_t b1;
  std::size_t b2;
};

const pt::idx_t EXPECTED_HAIRPIN_COUNT {1000}; // Expected number of hairpins
const boundary NULL_BOUNDARY{.b1 = pc::INVALID_IDX, .b2 = pc::INVALID_IDX};

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
} // namespace povu::algorithms
