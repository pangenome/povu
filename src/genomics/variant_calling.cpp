#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ctime>
#include <deque>
#include <format>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <stack>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "../common/typedefs.hpp"
//#include "../core/constants.hpp"
//#include "../core/core.hpp"
//#include "../common/utils.hpp"
#include "../graph/bidirected.hpp"
// #include "../graph/digraph.hpp"
//#include "../graph/tree.hpp"
//#include "../pvst/pvst.hpp"
#include "../io/io.hpp"


namespace genomics {
using namespace common::typedefs;




/**
 * Associate each path with a set of haplotypes
 *
 */
std::vector<std::vector<std::set<std::size_t>>>
find_path_haplotypes(
  const std::vector<std::vector<std::vector<bidirected::side_n_id_t>>>& all_paths,
  const bidirected::VariationGraph& bd_vg
  ) {
  std::string fn_name =  "[povu::genomics::find_path_haplotypes]";

  std::vector<std::vector<std::set<std::size_t>>> haplotypes_per_path;

  for (std::size_t fl_idx{}; fl_idx < all_paths.size(); ++fl_idx) {
    const std::vector<std::vector<bidirected::side_n_id_t>>& c_flubble_paths = all_paths[fl_idx];

    std::vector<std::set<std::size_t>> h;

    for (const std::vector<bidirected::side_n_id_t>& c_flubble_path: c_flubble_paths) {

      // for each path in the flubble
      // find the haplotype
      // add the haplotype to the path
      std::set<std::size_t> curr_haplotypes;
      std::set<std::size_t> intersect;
      std::set<std::size_t> temp;
      std::vector<std::vector<std::size_t>> haplotypes_per_path;
      for (auto [side, v_id]: c_flubble_path) {
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
 *
 *
 */
void call_variants(const std::vector<size_t_pair>& canonical_flubbles,
                   const bidirected::VariationGraph& bd_vg,
                   const core::config& app_config) {
  std::string fn_name = "[povu::genomics::call_variants]";
  if (app_config.verbosity() > 3) { std::cerr << fn_name << "\n"; }

  //output_format of = output_format::VCF;
  //std::vector<std::pair<std::size_t, std::size_t>> canonical_flubbles = extract_canonical_flubbles(pvst_, app_config);

  if (app_config.verbosity() > 3) {
    std::cerr << fn_name << " Found " << canonical_flubbles.size() << " canonical flubbles\n";
  }

  // walk paths in the digraph
  // while looping over canonical_flubbles
  std::vector<std::vector<std::vector<bidirected::side_n_id_t>>> all_paths;

  //std::cout << fn_name << " Extracting paths for flubbles:\n";
  // extract flubble paths
  for (std::size_t i{} ; i < canonical_flubbles.size(); ++i) {
    if (false) {
      std::cout << fn_name << " flubble: " << i
                << " start: " << canonical_flubbles[i].first
                << " stop: " << canonical_flubbles[i].second
                << "\n";
    }

    std::vector<std::vector<bidirected::side_n_id_t>> paths =
      bd_vg.get_paths(canonical_flubbles[i].first, canonical_flubbles[i].second);

    all_paths.push_back(paths);
  }

  //std::cout << fn_name << " Finding haplotypes\n";
  std::vector<std::vector<std::set<std::size_t>>> haplotypes_per_path =
    find_path_haplotypes(all_paths, bd_vg);

  std::map<std::size_t, std::vector<vcf::vcf_record>> vcf_records =
    vcf::gen_vcf_records(bd_vg, haplotypes_per_path, all_paths, app_config);

  vcf::write_vcfs(vcf_records, bd_vg, app_config);

}

} // namespace genomics
