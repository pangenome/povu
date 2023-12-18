#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iostream>
#include <stack>
#include <queue>
#include <string>
#include <map>
#include <set>
#include <unordered_set>
#include <vector>
#include <utility>
#include <ctime>
#include <iomanip>
#include <deque>
#include <format>
#include <sstream>
#include <fstream>

#include "./genomics.hpp"
#include "../pvst/pvst.hpp"
#include "../core/core.hpp"
#include "../graph/tree.hpp"
#include "../graph/digraph.hpp"
#include "../core/constants.hpp"
#include "../graph/bidirected.hpp"
#include "../core/utils.hpp"


namespace genomics {


/**
 * @brief get the edge indexes of a vertex on a given side
 */
const std::set<std::size_t>& get_edges_for_side(bidirected::Vertex const& v, bidirected::VertexEnd side) {
  return side == bidirected::VertexEnd::l ? v.get_edges_l() : v.get_edges_r();
}

/**
 * @brief compute the opposite side of an vertex
 * @return the opposite side of the given vertex
 */
bidirected::VertexEnd complement_side(bidirected::VertexEnd side) {
  return side == bidirected::VertexEnd::l ? bidirected::VertexEnd::r : bidirected::VertexEnd::l;
};


/**
 * @brief find the side of a vertex that is incident with one edge, if both are incident with one edge, return the right side
 *
 * @param bi_vg the bidirected variation graph
 * @param v_id the vertex id
 * @return the side that is incident with one edge, if both are incident with one edge, return the right side
 */
bidirected::VertexEnd compute_single_edge_side(const bidirected::VariationGraph& bi_vg, id_t v_id) {
  bidirected::Vertex const& v = bi_vg.get_vertex(v_id);
  return v.get_edges_r().size() == 1 ? bidirected::VertexEnd::r : bidirected::VertexEnd::l;
}

// TODO: this can be done while constructing the PVST
std::vector<std::pair<std::size_t, std::size_t>>
extract_canonical_flubbles(const tree::Tree& pvst_) {

  std::string fn_name{"[povu::genomics::extract_canonical_flubbles]"};

  if (true) {
	std::cerr << fn_name << std::endl;
  }
  
  tree::Vertex const& root = pvst_.get_root();

  // subtree set

  std::vector<std::pair<std::size_t, std::size_t>> canonical_flubbles;
  std::size_t current_vertex;

  // create a new queue
  std::queue<std::size_t> q;
  q.push(root.get_id());

  std::size_t max_child{core::constants::SIZE_T_MIN};
  bool is_canonical_subtree = true;
  while (!q.empty()) {
	current_vertex = q.front();
	q.pop();

	//std::cout << "current vertex: " << current_vertex << std::endl;

	std::set<std::size_t> const& children = pvst_.get_children(current_vertex);
	std::size_t parent_eq_class = pvst_.get_class(current_vertex);

	if (children.empty()) { continue; }

	is_canonical_subtree = true;
	bool end_found{false};

	max_child = core::constants::SIZE_T_MIN;

	for (std::size_t child: children) {
	  std::set<std::size_t> const& c = pvst_.get_children(child);
	  std::size_t child_eq_class = pvst_.get_class(child);

	  if (!end_found && child > max_child && parent_eq_class == child_eq_class) {
		max_child = child;
		end_found = true;
	  }

	  // all the children are leaves
	  if (!c.empty()) {
		q.push(child);
		is_canonical_subtree = false;
	  }
	}

	if (is_canonical_subtree) {
	  std::cerr << "canonical subtree found: " << current_vertex << " " << max_child
				<< " " << pvst_.get_meta(current_vertex) << " " << pvst_.get_meta(max_child) << std::endl;
	  if (pvst_.get_meta(current_vertex) < pvst_.get_meta(max_child)) {
		canonical_flubbles.push_back(
		  std::make_pair(pvst_.get_meta(current_vertex) - 1, pvst_.get_meta(max_child) - 1 )
		  );		
	  }
	  else {
		std::cerr << "sus\n";
	  }

	}
  }

  // loop over canonical subtrees and print
  /*
  std::cout << "Canonical subtrees:" << std::endl;
  for (auto subtree: canonical_flubbles) {
	std::cout << subtree.first << " " << subtree.second << std::endl;
  }
  */

  // use canonical subtrees to make a vector of pairs of flubble starts and stops
  //std::vector<std::pair<std::size_t, std::size_t>> canonical_flubbles;
  return canonical_flubbles;
}



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
void call_variants(const tree::Tree& pvst_,
				   const bidirected::VariationGraph& bd_vg,
				   const core::config& app_config) {
  std::string fn_name = "[povu::genomics::call_variants]";
  if (app_config.verbosity() > 3) { std::cerr << fn_name << "\n"; }

  //output_format of = output_format::VCF;
  std::vector<std::pair<std::size_t, std::size_t>> canonical_flubbles =
	extract_canonical_flubbles(pvst_);


  if (true) {
	std::cerr << fn_name << " Found " << canonical_flubbles.size() << " canonical flubbles\n";
  }
  //

  // walk paths in the digraph
  // while looping over canonical_flubbles


  std::vector<std::vector<std::vector<bidirected::side_n_id_t>>> all_paths;

  //std::cout << fn_name << " Extracting paths for flubbles:\n";
  // extract flubble paths
  for (std::size_t i{} ; i < canonical_flubbles.size(); ++i) {
	if (true) {
	  std::cout << fn_name << " flubble: " << i
				<< " start: " << canonical_flubbles[i].first
				<< " stop: " << canonical_flubbles[i].second
				<< "\n";
	}

	std::vector<std::vector<bidirected::side_n_id_t>>
	  paths = bd_vg.get_paths(canonical_flubbles[i].first,
							   canonical_flubbles[i].second);

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
