#include <cstddef>
#include <format>
#include <iostream>
#include <string>
#include <sys/types.h>
#include <unistd.h>
#include <vector>

#include "./cycle_equiv.hpp"
#include "../common/typedefs.hpp"



namespace algorithms {
using namespace graph_types;
using namespace common::typedefs;
using common::constants::INVALID_ID;


struct class_props {
  std::vector<std::size_t> g_id_; //
  std::vector<std::size_t> classes_; //
  std::vector<std::size_t> prev_; //
  std::vector<std::size_t> next_; //

  void print_data() {
    std::cout << "id: ";
      for (auto i : g_id_) {
        std::cout << i << " ";
      }
      std::cout << std::endl;

    std::cout << "c: ";
      for (auto i : classes_) {
        std::cout << i << " ";
      }
      std::cout << std::endl;


      // print prev
      std::cout << "prev: ";
      for (auto i : prev_) {
        std::cout << i << " ";
      }
      std::cout << std::endl;

            // print next
      std::cout << "next: ";
      for (auto i : next_) {
        std::cout << i << " ";
      }
      std::cout << std::endl;
  }
};

class_props stack_props(std::vector<id_n_cls> v) {
  std::string fn_name = std::format("[povu::genomics::graph_operations::{}]", __func__);

  std::map<std::size_t, std::size_t> eq_class_count;
  std::vector<std::size_t> eq_classes;
  std::vector<std::size_t> ids;

  for (std::size_t i{} ; i < v.size(); ++i) {
    if (eq_class_count.find(v[i].cls) == eq_class_count.end()) {
      eq_class_count[v[i].cls] = 1;
    }
    else {
      eq_class_count[v[i].cls] += 1;
    }
    ids.push_back(v[i].id);
    eq_classes.push_back(eq_class_count[v[i].cls]);
  }

  // a vector which contains the index within v of the next element with the same value
  std::map<std::size_t, std::size_t> last_seen;
  std::vector<std::size_t> nexter(v.size());
  for (std::size_t i{ v.size() } ; i-- > 0;) { // dec after checking
    if (last_seen.find(v[i].cls) == last_seen.end()) {
      nexter[i] = i;
    }
    else {
      nexter[i] = last_seen[v[i].cls];
    }
    last_seen[v[i].cls] = i;
  }

  last_seen.clear();
  std::vector<std::size_t> laster(v.size() );
  for (std::size_t i{} ; i < v.size(); ++i) {
    if (last_seen.find(v[i].cls) == last_seen.end()) {
      laster[i] = i;
    }
    else {
      laster[i] = last_seen[v[i].cls];
    }
    last_seen[v[i].cls] = i;
  }

  return {ids, eq_classes, laster, nexter};
}

std::vector<id_n_cls> compute_eq_class_stack(spanning_tree::Tree& t) {
  std::string fn_name = std::format("[povu::algorithms::{}]", __func__);
  std::vector<id_n_cls> s;

  std::set<std::size_t> seen;
  std::size_t g_v_id {};
  std::size_t curr_class {common::constants::UNDEFINED_SIZE_T};

  for (std::size_t v{t.size() - 1}; v < INVALID_ID; --v) {
    curr_class = common::constants::UNDEFINED_SIZE_T;

    try {
      g_v_id = std::stoull(t.get_vertex(v).name());
    }
    catch (const std::exception& e) {
      continue;
    }

    // g_v_id = std::stoull(t.get_vertex(v).name());
    if (t.is_root(v) || seen.find(g_v_id) != seen.end()) { continue; }
    seen.insert(g_v_id);

    const spanning_tree::Edge& e =  t.get_parent_edge(v);

    if (e.get_color() == color::black) {
      curr_class = e.get_class();
    }
    else {
      std::set<std::size_t> obe_idxs = t.get_obe_idxs(v);
      std::set<std::size_t> ibe_idxs = t.get_ibe_idxs(v);
      std::vector<std::size_t> be_idxs;

      be_idxs.reserve(obe_idxs.size() + ibe_idxs.size());
      be_idxs.insert(be_idxs.end(), obe_idxs.begin(), obe_idxs.end());
      be_idxs.insert(be_idxs.end(), ibe_idxs.begin(), ibe_idxs.end());

      for (auto be_idx : be_idxs) {
        if (t.get_backedge(be_idx).get_color() == color::black) {
          curr_class = t.get_backedge(be_idx).get_class();
          break;
        }
      }
    }

    if (curr_class != common::constants::UNDEFINED_SIZE_T) {
      s.push_back({g_v_id, curr_class});
    }
    else {
      throw std::runtime_error("No class found for vertex: " + t.get_vertex(v).name());
    }
  }

  if (false) {
    // print the contnets of s
    for (auto i : s) { std::cout << i.id << " "; }
    std::cout << std::endl;

    for (auto i : s) { std::cout << i.cls << " "; }
    std::cout << std::endl;
  }

  return s;
}

std::vector<size_t_pair> canonical_seses(const class_props& p) {
  std::vector<size_t_pair> v;
  const std::vector<std::size_t>& next { p.next_ };


  return v;
}

void find_seses(spanning_tree::Tree& t) {
  std::string fn_name = std::format("[povu::algorithms::{}]", __func__);
  std::vector<id_n_cls> s { compute_eq_class_stack(t) };
  class_props p { stack_props(s) };
  canonical_seses(p);
  //p.print_data();
}

} // namespace algorithms
