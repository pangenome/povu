#include <cstddef>
#include <iostream>
#include <limits>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <stack>
#include <utility>
#include <vector>



#include "./flubble_tree.hpp"
#include "../common/types.hpp"
#include "../common/utils.hpp"


namespace povu::graph::flubble_tree {

using namespace povu::graph_types;
namespace pc = povu::constants;
namespace pt = povu::types;
namespace pvtr = povu::tree;


// orientation, id, class
struct oic {
  pgt::or_t orientation;
  pt::id_t id;
  std::size_t cls;
};

inline flubble forwardise(pt::id_t start_id, pgt::or_t start_or, pt::id_t end_id, pgt::or_t end_or) {
  //flubble fl;
  if (start_or == pgt::or_t::reverse && end_or == pgt::or_t::reverse) {
    return {pgt::id_n_orientation_t{end_id, pgt::or_t::forward},
          pgt::id_n_orientation_t{start_id, pgt::or_t::forward}};
  }
  else {
    return {pgt::id_n_orientation_t{start_id, start_or},
          pgt::id_n_orientation_t{end_id, end_or}};
  }
}

/**
 * @brief Compute the equivalence class stack
 *
 * @param t The spanning tree
 * @return The equivalence class stack
 */
std::vector<oic> compute_eq_class_stack(pst::Tree& t) {
  std::string fn_name = std::format("[povu::algorithms::flubble_tree::{}]", __func__);

  // TODO: rename s to eq class stack
  std::vector<std::list<oic> *> s (t.size(), nullptr);
  std::set<std::size_t> seen;
  std::size_t g_v_id {}; // the sequence id of the vertex as in the GFA
  std::size_t curr_class { pc::UNDEFINED_SIZE_T };

  // "inherit" from child or create new entry in s
  auto add_to_stack = [&](std::size_t v) {
    const pst::Vertex& curr_vtx = t.get_vertex(v);

    if (curr_vtx.is_leaf()) {
      s[v] = new std::list<oic>{};
    }
    else if (curr_vtx.get_children().size() > 1) { // is a branching point
      s[v] = new std::list<oic>{};

      for (auto c : t.get_children(v)) {
        auto it = (curr_vtx.hi() == t.get_vertex(c).hi()) ? s[v]->end() : s[v]->begin();
        s[v]->splice(it, *s[c]);
        delete s[c];
        s[c] = nullptr;
      }
    }
    else { // linear and has one child
      std::size_t child_v_idx = *t.get_children(v).begin();

      s[v] = s[child_v_idx];
      s[child_v_idx] = nullptr;
    }
  };

  for (std::size_t v{t.size() - 1}; v < pc::INVALID_ID; --v) {
    curr_class = pc::UNDEFINED_SIZE_T;

    add_to_stack(v);

    pgt::v_type curr_vtx_type = t.get_vertex(v).type();

    if (curr_vtx_type != pgt::v_type::dummy) {
      g_v_id = std::stoull(t.get_vertex(v).name());
    }

    if (t.is_root(v) || seen.find(g_v_id) != seen.end()) { continue; }
    seen.insert(g_v_id);

    pgt::or_t curr_or;

    const pst::Edge& e =  t.get_parent_edge(v);
    if (e.get_color() == color::black) {
      curr_class = e.get_class();
      curr_or = curr_vtx_type == pgt::v_type::r ? pgt::or_t::forward : pgt::or_t::reverse;
    }
    else {
      bool found_black_edge { false };

      const std::set<std::size_t>& obe_idxs = t.get_obe_idxs(v);
      const std::set<std::size_t>& ibe_idxs = t.get_ibe_idxs(v);

      for (std::size_t be_idx : obe_idxs) {
        if (t.get_backedge(be_idx).get_color() == color::black) {
          curr_or = curr_vtx_type == pgt::v_type::r ? pgt::or_t::forward : pgt::or_t::reverse;
          found_black_edge = true;
          break;
        }
      }

      if (!found_black_edge) {
        for (std::size_t be_idx : ibe_idxs) {
          if (t.get_backedge(be_idx).get_color() == color::black) {
            curr_or = curr_vtx_type == pgt::v_type::r ? pgt::or_t::forward : pgt::or_t::reverse;
            found_black_edge = true;
            break;
          }
        }
      }
    }

    if (curr_class != pc::UNDEFINED_SIZE_T) {
      s[v]->push_front({curr_or, g_v_id, curr_class});
    }
    else {
      throw std::runtime_error(std::format("{} No class found for vertex: {}", fn_name, t.get_vertex(v).name()));
    }
  }

  std::vector<oic> stack_(s[t.get_root_idx()]->begin(), s[t.get_root_idx()]->end());
  delete s[t.get_root_idx()];

  return stack_;
}

/**
  * @brief Enumerate the flubbles
 */
pvtr::Tree<flubble> construct_flubble_tree(const std::vector<oic>& stack_, const std::vector<std::size_t>& next_seen) {
  std::string fn_name = std::format("[povu::algorithms::flubble_tree::{}]", __func__);

  pvtr::Tree<flubble> ft;

  struct ci {
    std::size_t cl; // class
    std::size_t idx; // index in stack_
  };
  std::stack<ci> s;
  std::unordered_set<std::size_t> in_s; // classes in s

  std::size_t prt_v { ft.root_idx() }; // parent vertex


  for (std::size_t i {}; i < stack_.size(); ++i) {
    auto [or_curr, id_curr, cl_curr] = stack_[i];

    if (in_s.count(cl_curr)) {

      while (!s.empty() && s.top().cl != cl_curr) {
        s.pop();
        in_s.erase(cl_curr);
      }

      s.pop();
      if (prt_v != ft.root_idx()) { prt_v = ft.get_parent_idx(prt_v); }
    } else {

    }

    if (i + 1 < next_seen[i]) {
      auto [or_nxt, id_nxt, _] = stack_[next_seen[i]];
      flubble fl = forwardise(id_curr, or_curr, id_nxt, or_nxt);
      pvtr::Vertex<flubble> v(id_curr, fl);
      std::size_t v_idx = ft.add_vertex(v);
      ft.add_edge(prt_v, v_idx);
      prt_v = v_idx;
    }

    s.push({cl_curr, i});
    in_s.insert(cl_curr);
  }

  return ft;
}

void compute_eq_class_metadata(const std::vector<oic> &stack_, std::vector<std::size_t>& next_seen) {

  std::unordered_map<std::size_t, std::size_t> next_seen_map;
  next_seen_map.reserve(stack_.size());

  for (std::size_t i {stack_.size()};  i-- > 0; ) {
    auto [or_curr, id_curr, cl_curr] = stack_[i];
    std::size_t next_idx = (next_seen_map.contains(cl_curr)) ? next_seen_map[cl_curr] : i;
    next_seen[i] = next_idx;
    next_seen_map[cl_curr] = i;
  }
}

/**
  * @brief Enumerate the flubbles
 */
std::vector<flubble> find_flubbles(const std::vector<oic>& stack_) {
  std::string fn_name = std::format("[povu::algorithms::flubble_tree::{}]", __func__);

  std::vector<flubble> flubbles;
  // class to last seen idx
  std::unordered_map<std::size_t, std::size_t> last_seen_map;
  // the index of the last time the class at index i in stack_ was seen
  std::vector<std::size_t> last_seen;
  last_seen.reserve(stack_.size());

  for (std::size_t i{}; i < stack_.size(); ++i) {
    auto [or_curr, id_curr, cl_curr] = stack_[i];
    std::size_t prev_idx = (last_seen_map.contains(cl_curr)) ? last_seen_map[cl_curr] : i;
    last_seen.push_back(prev_idx);

    if (prev_idx + 1 < i) {
      auto [or_prev, id_prev, cl_prev] = stack_[prev_idx];
      flubble fl = forwardise(id_prev, or_prev, id_curr, or_curr);
      flubbles.push_back(fl);
    }

    last_seen_map[cl_curr] = i;
  }

  return flubbles;
}

pvtr::Tree<flubble> st_to_ft(pst::Tree& t) {
  std::string fn_name = std::format("[povu::algorithms::{}]", __func__);

  std::vector<oic> s{compute_eq_class_stack(t)};

  std::vector<std::size_t> next_seen (s.size(), pc::INVALID_IDX);
  compute_eq_class_metadata(s, next_seen);


  pvtr::Tree<flubble> ft = construct_flubble_tree(s, next_seen);

  return ft;
}

std::vector<flubble> enumerate(pst::Tree& t) {
  std::string fn_name = std::format("[povu::algorithms::{}]", __func__);

  std::chrono::duration<double> timeRefRead;
  auto t0 = pt::Time::now();

  std::vector<oic> s { compute_eq_class_stack(t) };

  bool time {true};

  if (time) {
    timeRefRead = pt::Time::now() - t0;
    povu::utils::report_time(std::cerr, fn_name, "compute equiv class stack", timeRefRead);
    t0 = pt::Time::now();
  }

  std::vector<flubble> flubbles = find_flubbles(s);

  if (time) {
    timeRefRead = pt::Time::now() - t0;
    povu::utils::report_time(std::cerr, fn_name, "enumerate flubbles", timeRefRead);
    t0 = pt::Time::now();
  }

  return flubbles;
}

}
