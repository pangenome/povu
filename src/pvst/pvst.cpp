#include <algorithm>
#include <bit>
#include <cmath>
#include <cstddef>
#include <format>
#include <iostream>
#include <map>
#include <string>
#include <sys/types.h>
#include <vector>
#include <list>
#include <fstream>

#include "./pvst.hpp"
#include "../graph/tree.hpp"
#include "../core/core.hpp"
#include "../core/utils.hpp"

using namespace core::constants;

namespace pvst {


/**
 * @brief compute the pvst of a given vector of eq_n_id_t values
 *
 * @param v a vector of pairs where each pair is (index, eq class)
 * @param app_config the application configuration
 *
 * @return a PVST tree
 */
tree::Tree compute_pvst(std::vector<core::eq_n_id_t> v, const core::config& app_config) {
  std::string fn_name = "[povu::pvst::compute_pvst]";
  if (app_config.verbosity() > 3) { std::cerr << fn_name << "\n"; }

  tree::Tree t = tree::Tree();

  // lambda to map back an index from the biedged and dummies vertex idx to the
  // bidirected idx
  auto bidirected_idx = [](std::size_t x) -> std::size_t {
    --x; // we added 1 because we added a dummy start node before bi-edging
    return x % 2 == 0 ? ((x + 2) / 2) - 1 : ((x + 1) / 2) - 1;
  };

  // a map of equivalence classes to a vector of their positions in the edge stack
  std::map<size_t, std::vector<size_t>> pos_map;
  std::map<size_t, bool> nesting;
  for (std::size_t i{}; i < v.size(); ++i) {
    pos_map[v[i].eq_class].push_back(i);
    nesting[v[i].eq_class] = false;
  }

  // print pos_map
  if (false) {
    std::cout << "pos_map\n";
    for (auto i: pos_map) {
      std::cout << i.first << ": ";
      for (auto j: i.second) { std::cout << j << ", "; }
      std::cout << std::endl;
    }
  }

  std::string bd_idx_str{};

  std::size_t
    current_class{}, //
    counter{}, // a counter to keep track of the number of vertices in the tree
    bd_idx{}, // bidirected index
    be_idx{},  // biedged index
    i{}, // a position in the edge stack (v)
    p_id{INVALID_ID}, // parent id of the current vertex in the tree
    p_class // parent class
    ;

  tree::Vertex p_v; // parent vertex
  // tree::Vertex prt_v;

  auto is_root =[&](id_t id) -> bool {
    return p_id == INVALID_ID || t.get_root().get_id() == id;
  };

  while (i < v.size()) {
    be_idx = v[i].v_id;
    bd_idx = bidirected_idx(be_idx);
    current_class = v[i].eq_class;
    bd_idx_str = std::to_string(bd_idx+1);

    if (false) {
       std::cerr << std::format("i: {} counter {} p_id: {} be idx: {} bd idx: {} current class: {}\n",
                                i, counter, p_id, be_idx, bd_idx, current_class);
    }

    // TODO: is it necessary to check if i == 0? or the state of the tree?

    if (t.empty()) { // add root
      t.add_vertex(INVALID_ID, counter, current_class, bd_idx_str);
    }
    else { // add other vertices
      t.add_vertex(p_id, counter, current_class, bd_idx_str);
    }

    // ---------------------------
    // update/determine the parent
    // ---------------------------

    std::vector<std::size_t> const& ps = pos_map[current_class];

    // if this eq class shows up more than once
    if (ps.size() > 1) {
      // in this case the eq class nests something or is in a bubble chain
      std::list<std::size_t> positions;
      std::copy(ps.begin(), ps.end(), std::back_inserter(positions));

      for (auto it = positions.begin(); it != positions.end(); ++it) {

        // detect & handle nesting and bubble chains
        if (*it == i && std::next(it) != positions.end() && it != positions.begin() && ((*std::next(it) - *it ) > 1) ) {
          // std::cout << "nesting\n";
          tree::Vertex& prt_v = t.get_vertex_mut(p_id);

          if (!prt_v.is_dummy()) {
            if (is_root(p_id)) {
              // std::cerr << "prt is root\n";
              std::string m = prt_v.get_meta();
              ++counter;
              t.add_vertex(INVALID_ID, counter, prt_v.get_class(), m, true);
              t.set_root(counter);
              // std::cerr << "p id " << p_id << " ctr " << counter << std::endl;
              t.get_vertex_mut(p_id).set_parent(counter);
              prt_v.set_parent(counter);
              // std::cout << prt_v.get_parent() << " p " << t.get_vertex(p_id).get_parent() << std::endl;

              // std::cerr << "children count: " << t.get_vertex_mut(p_id).get_children().size() << std::endl;

              auto children = t.get_vertex_mut(p_id).get_children();

              for (auto c:  children) {
                t.get_vertex_mut(c).set_parent(counter);
                t.get_vertex_mut(counter).add_child(c);
              }

              for (auto c:  children) {
                t.get_vertex_mut(p_id).remove_child(c);
              }

              t.get_vertex_mut(counter).add_child(p_id);
              t.add_vertex(p_id, counter, current_class, bd_idx_str, true);

              p_class = prt_v.get_class();
            }
            else {
              // std::cerr << "prt not root\n";
              id_t grand_p_id = t.get_parent(prt_v.get_id());
              tree::Vertex& grand_prt_v = t.get_vertex_mut(grand_p_id);

              t.get_vertex_mut(grand_p_id).remove_child(p_id);
              std::string m = grand_prt_v.get_meta();
              ++counter;
              t.add_vertex(grand_p_id, counter, grand_prt_v.get_class(), m, true);
              t.get_vertex_mut(counter).add_child(p_id);

              p_class = grand_prt_v.get_class();
            }

            p_id = counter;
          }
          else {
            if (!is_root(prt_v.get_id())) {
              p_id = prt_v.get_parent();
            }
          }

          ++counter;
          t.add_vertex(p_id, counter, current_class, bd_idx_str, true);
          nesting[current_class] = false;
        }
      }
    }

    bool foo{false}; // when parent and child are in the same eq class

    if ((!nesting[current_class] && ps.size() > 1) || (!nesting[current_class] && i+1 < v.size() && current_class != v[i+1].eq_class)) {
      // std::cout << "down\n";
      nesting[current_class] = true;
      if (p_class == current_class) { foo = true; }
      p_class = current_class;
      p_id = counter;
    }

    if (ps.back() == i || (!t.get_vertex(counter).is_dummy() && foo && i+1 < v.size() && current_class != v[i+1].eq_class)) {
      // std::cout << "flip\n";
      nesting[current_class] = false;
    }

    if (!nesting[p_class] && !is_root(p_id)) {
      // std::cout << "up\n";
      // this is the last time we are seeing this class; go up // get the parent's parent
      p_id = t.get_parent(p_id);
      p_v = t.get_vertex(p_id);
      p_class = p_v.get_class();
    }

    ++i;
    ++counter;
  }

  return t;
}

void to_text(tree::Tree const& pvst, std::filesystem::path const& output_path) {
  std::ofstream f(output_path);

  if (!f.is_open()) {
    std::cerr << "ERROR: could not open file " << output_path << "\n";
    std::exit(1);
  }

  for (id_t i{0}; i < pvst.size(); ++i) {
    const tree::Vertex& v = pvst.get_vertex(i);

    std::set<std::size_t> c = v.get_children();
    std::set<std::string> c_ {};
    std::transform(c.begin(), c.end(), std::inserter(c_, c_.begin()),
                   [](int value) { return std::to_string(value); });

    f << v.get_id()  << "\t";
    f << (utils::concat_with(c_,  ',')) << "\t";
    f << (v.is_dummy() ? "D" : "T") << "\t";
    f << v.get_class() << "\t";
    f << v.get_meta() << "\n";
  }

  f.close();
}
} // namespace pvst
