#include <algorithm>
#include <cstddef>
#include <ctime>
#include <format>
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "../graph/bidirected.hpp"
#include "../io/io.hpp"


namespace genomics {
typedef std::vector<bidirected::id_n_orientation_t> walk; // a walk is a sequence of vertices also a path

/**
 * Associate each path with a set of haplotypes
 *
 */
std::vector<std::vector<std::set<std::size_t>>>
find_path_haplotypes(const std::vector<std::vector<walk>>& all_paths, const bidirected::VariationGraph& bd_vg) {
  std::string fn_name { std::format("[povu::genomics::{}]",  __func__) };

  std::vector<std::vector<std::set<std::size_t>>> haplotypes_per_path;

  for (std::size_t fl_idx{}; fl_idx < all_paths.size(); ++fl_idx) {
    const std::vector<walk>& c_flubble_paths = all_paths[fl_idx];

    std::vector<std::set<std::size_t>> h;

    for (const walk& c_flubble_path: c_flubble_paths) {

      // for each path in the flubble
      // find the haplotype
      // add the haplotype to the path
      std::set<std::size_t> curr_haplotypes;
      std::set<std::size_t> intersect;
      std::set<std::size_t> temp;
      std::vector<std::vector<std::size_t>> haplotypes_per_path;
      for (auto [v_id, _]: c_flubble_path) {
        // find the haplotype
        const std::vector<bidirected::PathInfo>& path_info = bd_vg.get_vertex(v_id).get_paths();


        for (auto [path_id, _]: path_info) {
          curr_haplotypes.insert(path_id);
        }

        if (intersect.empty()) {
          intersect = curr_haplotypes;
        }

        std::set_intersection(curr_haplotypes.begin(), curr_haplotypes.end(),
                              intersect.begin(), intersect.end(),
                              std::inserter(temp, temp.begin()));
        intersect = temp;
        temp.clear();
      }

      h.push_back(intersect);
    }

    haplotypes_per_path.push_back(h);
  }

  return haplotypes_per_path;
}

/**
 * @brief Get the side of a start or end vertex that is incident to the SESE
 *
 * the alternative side faces the SESE
 * if the incident side is the left side then the vertex is in forward orientation and vice versa
 * Being an SESE we expect that all vertices connected to one side are in the SESE
 * and all vertices connected to the other side are not in the SESE
 *
 * @param in_sese The set of vertices in the SESE
 * @param start_id The id of the vertex to get the orientation for
 */
bidirected::VertexEnd get_boundary_incidence(const bidirected::VariationGraph& g,
                                 const std::set<std::size_t>& in_sese,
                                 std::size_t v_id,
                                 std::size_t alt_id) {
  std::string fn_name = std::format("[povu::genomics::]", __func__);

  std::size_t v_idx = g.id_to_idx(v_id);
  const bidirected::Vertex& v = g.get_vertex(v_idx);

  auto foo = [&](std::size_t e_idx) ->bool {
    auto [side, alt_v_idx] = g.get_edge(e_idx).get_other_vertex(v_idx);
    return in_sese.count(g.idx_to_id(alt_v_idx)) || g.idx_to_id(alt_v_idx) == alt_id ;
  };

  bool allLeftInSet = std::any_of(v.get_edges_l().begin(), v.get_edges_l().end(), foo);

  bool allRightNotInSet = allLeftInSet ? !std::none_of(v.get_edges_r().begin(), v.get_edges_r().end(), foo)
                                       : std::any_of(v.get_edges_r().begin(), v.get_edges_r().end(), foo);
  if (!(allLeftInSet ^ allRightNotInSet)) {
    throw std::runtime_error(std::format("{} {} {}", v_idx, allLeftInSet, allRightNotInSet));
  }

  bidirected::VertexEnd o = allLeftInSet ? bidirected::VertexEnd::l : bidirected::VertexEnd::r;

  return o;
}

std::pair<bidirected::id_n_orientation_t, bidirected::id_n_orientation_t>
foo(const bidirected::VariationGraph& g, const graph_types::canonical_sese& sese){
  std::string fn_name = std::format("[povu::genomics::{}]", __func__);
  if (false) { std::cerr << fn_name << "\n"; }

  auto [start_id, stop_id, in_sese] = sese;

  /*
    assume e is the incoming vertex end
   */
  auto v_end_to_orientation = [&](graph_types::VertexEnd e) -> bidirected::orientation_t {
    if (e == graph_types::VertexEnd::l ) {
      return bidirected::orientation_t::forward;
    }

    return bidirected::orientation_t::reverse;
  };

  /*
    Determine start side and orientation
   */
  bidirected::VertexEnd start_end;
  try {
    start_end = get_boundary_incidence(g, in_sese, start_id, stop_id);
  }
  catch (std::exception& e) {
    std::cerr << std::format("{} WARN: SESE ({} {}) start boundary: {}\n", fn_name, start_id, stop_id, e.what());
    return{};
  }

  bidirected::VertexEnd stop_end;
  try {
    stop_end = get_boundary_incidence(g, in_sese, stop_id, start_id);
  }
  catch (std::exception& e) {
    std::cerr << std::format("{} WARN: SESE ({} {}) stop boundary: {}\n", fn_name, start_id, stop_id, e.what());
    return{};
  }


  //std::cerr << start_orientation << " " << stop_orientation << "\n";

    // TODO: remove when sort is implemented
  // if SESE is flipped then we need to reverse the start and stop
  if (start_id > stop_id && start_end == bidirected::VertexEnd::l && stop_end == bidirected::VertexEnd::r) {
    std::swap(start_id, stop_id);
    std::swap(start_end, stop_end);
  }

  bidirected::orientation_t start_o = start_end == graph_types::VertexEnd::r ? bidirected::orientation_t::forward
                                                                            : bidirected::orientation_t::reverse;

  bidirected::orientation_t stop_o = stop_end == graph_types::VertexEnd::l ? bidirected::orientation_t::forward
                                                                          : bidirected::orientation_t::reverse;

  bidirected::id_n_orientation_t entry { start_id, start_o };
  bidirected::id_n_orientation_t exit { stop_id, stop_o };

  return { entry, exit };
}

/**
 *
 *
 */
void call_variants(const std::vector<graph_types::canonical_sese>& canonical_flubbles,
                   const bidirected::VariationGraph& bd_vg,
                   const core::config& app_config) {
  std::string fn_name { std::format("[povu::genomics::{}]" , __func__) };

  // walk paths in the digraph
  // while looping over canonical_flubbles
  std::vector<std::vector<walk>> all_paths;

  // extract flubble paths
  for (std::size_t i{}; i < canonical_flubbles.size(); ++i) {
    const graph_types::canonical_sese& f = canonical_flubbles[i];
    auto [entry, exit] = foo(bd_vg, f);
    std::vector<walk> paths = bd_vg.get_paths(entry, exit);
    all_paths.push_back(paths);
  }

  std::vector<std::vector<std::set<std::size_t>>> haplotypes_per_path =
    find_path_haplotypes(all_paths, bd_vg);


  std::map<std::size_t, std::vector<vcf::vcf_record>> vcf_records =
    vcf::gen_vcf_records(bd_vg, haplotypes_per_path, all_paths, app_config);

  vcf::write_vcfs(vcf_records, bd_vg, app_config);
}

} // namespace genomics
