#include <cstddef>
#include <ctime>
#include <format>
#include <iostream>
#include <map>
#include <optional>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "./genomics.hpp"
#include "../core/constants.hpp"
#include "../graph/tree.hpp"
#include "../common/common.hpp"


namespace graph_operations {
using namespace graph_types;

struct id_n_cls {
    std::size_t id;
    std::size_t cls;
};

  // implement << operator for id_n_cls
std::ostream& operator<<(std::ostream& os, const id_n_cls& id_cls) {
    os << "id: " << id_cls.id << " cls: " << id_cls.cls;
    return os;
}

// variant paths are the nodes whose children are only leaves
std::vector<std::size_t> get_variant_paths(const tree::MapTree<std::size_t, node>& pst) {
  // map of parent to children
  std::map<std::size_t, std::set<std::size_t>> variant_paths;
  std::set<std::size_t> leaf_parents;

  for (auto [v_idx, n] : pst) {
    if (n.children.empty() && n.parent != core::constants::UNDEFINED_SIZE_T) {
      variant_paths[n.parent].insert(v_idx);
      leaf_parents.insert(n.parent);
    }
  }

  // filter for leaf only parents
  for (std::size_t p : leaf_parents) {
    // have one child
    if (pst.at(p).children.size() == 1) {
      variant_paths.erase(p);
    }
    else {
      // if any child is not a leaf, remove the parent from variant paths
      for (std::size_t c : pst.at(p).children) {
        if (pst.at(c).children.size() > 0) {
          variant_paths.erase(p);
          break;
        }
      }
    }
  }

  // Extract keys and push them into the vector using transform
  std::vector<std::size_t> keys;
  std::transform(variant_paths.begin(),
                 variant_paths.end(),
                 std::back_inserter(keys),
                 [](const auto& p) { return p.first; });

  return keys;
}


/**
 * @brief print the parent child relationships (tree) in dot format as an undirected graph
 */
void pst_to_dot(const tree::MapTree<std::size_t, node>& pst) {
  std::cout << "graph G {\n"
            << "\trankdir = LR;\n"
            << "\tnode[shape = circle];\n"
            << "\tedge [arrowhead=vee];\n";

  // print a set of pairs with commas except for the last one

  auto vec_to_str = [](const std::set<std::pair<std::size_t, std::size_t>>& vec) {
    std::string str;
    for (auto [a, b] : vec) {
      str += std::format("({}, {}), ", a, b);
    }
    // delete the last ', '
    return str.substr(0, str.size() - 2);
      //return str;
  };

  for (auto [k, n] : pst) {
    std::cout << std::format("\t{} [label = \"{} {}\"];\n", k, k, vec_to_str(n.ids)) ;
  }

  for (auto [v_id, node] : pst) {
    for (auto child : node.children) {
      std::cout << std::format("\t{} -- {};\n", v_id, child);
    }
  }

  std::cout << "}" << std::endl;
}


/**
 * @brief not really a PST but a tree of canonical SESEs or nesting SESEs
 */
tree::MapTree<std::size_t, node>
compute_pst(const std::vector<id_n_cls> &v,
            std::map<std::size_t, std::size_t> eq_class_count,
            std::vector<std::size_t> nexter,
            std::vector<std::size_t> laster) {
  tree::MapTree<std::size_t, node> pst;
  std::size_t curr_class { core::constants::UNDEFINED_SIZE_T };

  // the last element with this value: key is the class, value is the index
  std::map<std::size_t, std::size_t> lasts;

  for (std::size_t i{} ; i < v.size(); ++i) {
    if (nexter[i] == i) {
      // this is the last element with this value last time we see this class
      // because it occurs once or it is the last element with this value
      if (eq_class_count[v[i].cls] > 1) {
        curr_class = pst[v[i].cls].parent;

        // if we don't have two consecutive classes at once
        if (laster[i] != i-1 || eq_class_count[v[i].cls] == 2) {
          pst[v[i].cls].ids.insert( {v[laster[i]].id, v[i].id} );
        }
      }
    }
    else {
      if (curr_class == core::constants::UNDEFINED_SIZE_T) {
        pst[v[i].cls] = {};
      }
      else if (curr_class != v[i].cls) {
        pst[curr_class].ids.insert( {v[lasts[curr_class]].id, v[nexter[lasts[curr_class]]].id} );

         pst[curr_class].children.insert(v[i].cls);
        pst[v[i].cls] = {};
        pst[v[i].cls].parent = curr_class;
      }

      curr_class = { v[i].cls };
    }

    if (eq_class_count[v[i].cls] > 1) { lasts[v[i].cls] = i; }
  }


  return pst;
}

tree::MapTree<std::size_t, node> find_seses(spanning_tree::Tree st, const std::vector<std::size_t>& edge_ids) {
  std::string fn_name = std::format("[povu::genomics::graph_operations::{}]", __func__);

  std::vector<id_n_cls> v;
  for (auto e_id : edge_ids) {
    auto [e_type, e_idx] = st.get_edge_idx(e_id);
    std::size_t cls = (e_type == spanning_tree::EdgeType::tree_edge) ?
      st.get_tree_edge(e_idx).get_class() : st.get_backedge(e_idx).get_class();
    v.push_back({e_id, cls});
  }

  std::map<std::size_t, std::size_t> eq_class_count;
  std::vector<std::size_t> eq_classes;

  for (std::size_t i{} ; i < v.size(); ++i) {
    if (eq_class_count.find(v[i].cls) == eq_class_count.end()) {
      eq_class_count[v[i].cls] = 1;
    }
    else {
      eq_class_count[v[i].cls] += 1;
    }
    eq_classes.push_back(eq_class_count[v[i].cls]);
  }

  // a vector which contains the index within v of the next element with the same value
  std::map<std::size_t, std::size_t> last_seen;
  std::vector<std::size_t> nexter(v.size() );
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

  if(false){
    std::cout << fn_name << std::endl;
    std::cout << "id: ";
    for (std::size_t i{} ; i < v.size(); ++i) {
      std::cout << v[i].id << " ";
    }
    std::cout << "\n";

    std::cout << "c: ";
    for (std::size_t i{} ; i < v.size(); ++i) {
      std::cout << v[i].cls << " ";
    }
    std::cout << "\n";

    std::cout << "eq_classes: ";
    for (std::size_t i{} ; i < eq_classes.size(); ++i) {
      std::cout << eq_classes[i] << " ";
    }
    std::cout << "\n";

    std::cout << "eq_class_count: ";
    for (std::size_t i{} ; i < v.size(); ++i) {
      std::cout << eq_class_count[v[i].cls] << " ";
    }
    std::cout << "\n";

    std::cout << "nexter: ";
    for (std::size_t i{} ; i < v.size(); ++i) {
      std::cout << nexter[i] << " ";
    }
    std::cout << "\n";
  }

  tree::MapTree<std::size_t, node> pst = compute_pst(v, eq_class_count, nexter, laster);

  if (true) {
    std::cout << std::format("{} PST\n", fn_name);
    //pst_to_dot(pst);
    std::cout << "\n";
  }

  return pst;
}

std::pair<std::size_t, std::size_t> extract_pair(spanning_tree::Tree st, std::size_t edge_id) {
  auto [start_e_type, start_e_idx] = st.get_edge_idx(edge_id);

  std::pair<std::size_t, std::size_t> k;

  if ( start_e_type == spanning_tree::EdgeType::tree_edge ) {
    const spanning_tree::Edge& start_edge = st.get_tree_edge(start_e_idx);
    k = st.get_vertex(start_edge.get_v1()).dfs_num() < st.get_vertex(start_edge.get_v2()).dfs_num() ?
      std::make_pair(start_edge.get_v1(), start_edge.get_v2()) : std::make_pair(start_edge.get_v2(), start_edge.get_v1());
  }
  else {
    const spanning_tree::BackEdge& start_edge = st.get_backedge(start_e_idx);
    k = {start_edge.get_src(), start_edge.get_tgt()};
  }

  return k;
}

std::vector<std::tuple<std::size_t, VertexType, std::size_t, VertexType>>
foo(spanning_tree::Tree st, const std::vector<std::size_t>& edge_ids, const core::config& app_config) {
  std::string fn_name = std::format("[povu::genomics::graph_operations::{}]", __func__);

  tree::MapTree<std::size_t, node> pst = find_seses(st, edge_ids);

  // each member of this vector is an SESE with variants
  std::vector<std::size_t> v_ = get_variant_paths(pst);

  std::vector<std::tuple<std::size_t, VertexType, std::size_t, VertexType>> seses;

  for (auto i : v_) {
    const std::set<std::pair<std::size_t, std::size_t>>& x = pst[i].ids;
    for (auto [start_id, end_id] : x) {

      std::pair<std::size_t, std::size_t> start_pair = extract_pair(st, start_id);
      std::pair<std::size_t, std::size_t> end_pair = extract_pair(st, end_id);

      auto [_, s] = start_pair;
      auto [e, __] = end_pair;

      st.get_vertex(s).type();
      st.get_vertex(e).type();

      std::size_t bd_s =  common_fns::to_bidirected_idx(s);
      std::size_t bd_e = common_fns::to_bidirected_idx(e);

      seses.push_back({bd_s, st.get_vertex(s).type(), bd_e, st.get_vertex(e).type()});

      //  std::cout << std::format("start: {}{} end : {}{}\n", s, bd_s, e, bd_e);
    }
  }

  return seses;
}

}
