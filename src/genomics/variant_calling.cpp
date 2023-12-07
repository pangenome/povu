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
 * @brief Given a bidirected::VariationGraph and a subpath, return the string sequence
 */
std::string path_to_seq(const bidirected::VariationGraph& bd_vg,
						const std::vector<bidirected::side_n_id_t>& p) {
  std::string seq = "";
  for (auto [side, id]: p) {
	bidirected::Vertex const& v = bd_vg.get_vertex(id);
	seq += side == bidirected::VertexEnd::l ? v.get_label() : v.get_rc_label();
  }

  return seq;
}

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

	for (auto child: children) {
	  std::set<std::size_t> const& c =  pvst_.get_children(child);
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
	  canonical_flubbles.push_back(
		std::make_pair(pvst_.get_meta(current_vertex) - 1, pvst_.get_meta(max_child) - 1 )
		);
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

// TODO: move to vcf.hpp
/**
 * @brief
 */
std::map<std::size_t, std::vector<vcf::vcf_record>>
gen_vcf_records(
  const bidirected::VariationGraph& bd_vg,
  const std::vector<std::vector<std::set<std::size_t>>>& haplotypes,
  const std::vector<std::vector<std::vector<bidirected::side_n_id_t>>>& all_paths,
  const core::config& app_config
  ) {

  // all our VCF files will be in this map
  std::map<std::size_t, std::vector<vcf::vcf_record>> vcf_records;

  std::vector<std::string> reference_paths = app_config.get_reference_paths();

  utils::TwoWayMap<std::size_t, std::string> paths_map; // path id to its name

  for (auto p : bd_vg.get_paths()) {
	paths_map.insert(p.id, p.name);
  }

  if (app_config.gen_undefined_vcf()) {
	// add the undefined reference for flubbles that don't have the a reference
	paths_map.insert(vcf::UNDEFINED_PATH_ID, vcf::UNDEFINED_PATH_LABEL);
	reference_paths.push_back(vcf::UNDEFINED_PATH_LABEL);
  }

  //std::size_t num_flubbles = haplotypes.size();
  std::size_t num_paths{};
  // for each flubble
  for (std::size_t fl_idx{}; fl_idx < all_paths.size(); ++fl_idx) {
	num_paths = all_paths[fl_idx].size();

	// does this flubble contain a reference path?
	bool has_ref{false};

	// while being strand aware
	// spell all the sequences in the flubble as a vector of strings
	// where each string is a path
	std::vector<std::string> path_seqs;
	for (std::size_t path_idx{}; path_idx < num_paths; ++path_idx) {
	  const std::vector<bidirected::side_n_id_t>& p = all_paths[fl_idx][path_idx];
	  path_seqs.push_back(path_to_seq(bd_vg, p));
	}

	// for each path in the flubble
	// if the path is not a reference path
	// then it is a variant
	// and we need to add it to the vcf
	// we look at the haplotypes that pass through the path
	for (std::size_t path_idx{}; path_idx < num_paths; ++path_idx) {

	  // get the haplotypes that pass through the path at path_idx
	  std::set<std::size_t> path_haplotypes = haplotypes[fl_idx][path_idx];

	  // not every path will be supported by at least one haplotype
	  if (path_haplotypes.empty()) {
		//std::cerr << "ERROR: path " << path_idx << " in flubble " << fl_idx << " has no haplotypes\n";
		continue;
	  }

	  // because we want a VCF for each reference path then
	  // for each reference path we check if it is in the set of haplotypes
	  // for that path
	  // if it is then remove it from the path seqs and the rest are variants
	  // we then populate the vcf records map which contains
	  // the vcf records for each reference path
	  for (const std::string& ref_name : reference_paths) {
		std::size_t ref_id = paths_map.get_key(ref_name);

		if (path_haplotypes.count(ref_id)) {

		  has_ref = true;

		  // the path contains the reference haplotype ref_id
		  vcf::vcf_record vcf_rec;

		  // TODO: this is a hack, a path may differ internally by a little
		  // use the intersection to determine this

		  // get the position of that path in the reference whose id is ref_id
		  auto [_, v_id] = all_paths[fl_idx][path_idx].front();
		  bidirected::Vertex const& v = bd_vg.get_vertex(v_id);

		  // TODO: can this loop be avoided?
		  for (const bidirected::PathInfo& vertex_path : v.get_paths()) {

			if (vertex_path.path_id == ref_id) {
			  vcf_rec.pos = vertex_path.step_index;
			}

			vcf_rec.ref = path_seqs[path_idx];
			vcf_rec.alt = utils::immutable_erase(path_seqs, path_idx);

			vcf_records[ref_id].push_back(vcf_rec);
		  }
		}
	  }
	}

	if (app_config.gen_undefined_vcf() && !has_ref) {
	  vcf::vcf_record vcf_rec;
	  vcf_rec.pos = vcf::UNDEFINED_PATH_POS;
	  vcf_rec.ref = "";
	  vcf_rec.alt = path_seqs;

	  vcf_records[vcf::UNDEFINED_PATH_ID].push_back(vcf_rec);
	}
  }

  return vcf_records;
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


  //std::cout << fn_name << " Found " << canonical_flubbles.size() << " canonical flubbles\n";

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

	std::vector<std::vector<bidirected::side_n_id_t>>
	  paths = bd_vg.get_paths(canonical_flubbles[i].first,
							   canonical_flubbles[i].second);

	all_paths.push_back(paths);
  }

  //std::cout << fn_name << " Finding haplotypes\n";
  std::vector<std::vector<std::set<std::size_t>>> haplotypes_per_path =
	find_path_haplotypes(all_paths, bd_vg);

  std::map<std::size_t, std::vector<vcf::vcf_record>> vcf_records =
	gen_vcf_records(bd_vg, haplotypes_per_path, all_paths, app_config);

  vcf::write_vcfs(vcf_records, bd_vg, app_config);

}

} // namespace genomics
