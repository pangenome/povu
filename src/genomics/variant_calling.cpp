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



namespace genomics {

typedef std::size_t id_t; //
typedef id_t vertex_id_t; //
typedef std::pair<bidirected::VertexEnd, id_t> side_n_id_t;
typedef std::vector<side_n_id_t> subpath_t;
typedef std::vector<subpath_t> subpaths_t;

const std::string UNDEFINED_PATH_LABEL = "undefined";
std::size_t UNDEFINED_PATH_ID = core::constants::SIZE_T_MAX;
std::size_t UNDEFINED_PATH_POS = core::constants::SIZE_T_MAX;
  
/**
 * @brief
 *
 * @param v the vector whose value is to be erased
 * @param idx the index to be erased
 * @return a vector with the value at the given index erased
 */
std::vector<std::string> immutable_erase(std::vector<std::string>& v, std::size_t idx) {
  std::vector<std::string> v_ = v;
  v_.erase(v_.begin()+idx, v_.begin()+idx);
  return v_;
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
	  canonical_flubbles.push_back(std::make_pair(current_vertex, max_child));
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
 * @brief Simplify paths by removing redundant nodes and instead representing it as a single node and strand through which the path passes
 *
 * @param all_paths a vector of paths
 * @return std::vector<subpaths_t> a vector of simplified paths
 */
subpaths_t simplify_paths(const subpaths_t& all_paths) {

  subpaths_t simplified_paths;
  simplified_paths.reserve(all_paths.size());

  for (std::size_t i{}; i < all_paths.size(); ++i) {
	const subpath_t& path = all_paths[i];

	if (path.size() % 2 != 0) {
	  std::cout << "Expected even number of paths, skipping path\n";
	  continue;
	}

	std::vector<side_n_id_t> simplified_path;
	simplified_path.reserve(path.size()/2);

	for (std::size_t j{}; j < path.size(); j+=2) {
	  if (path[j].second != path[j+1].second) {
		std::cout << "Expected same node id, skipping path\n";
	  }

	  if (path[j].first == path[j+1].first) {
		std::cout << "Expected different sides, skipping path\n";
	  }

	  simplified_path.push_back(std::make_pair(path[j].first, path[j].second));
	}

	simplified_paths.push_back(simplified_path);
  }

  return simplified_paths;
}

/**
 * @brief Get all paths between two nodes in the bidirected variation graph
 *
 * @param start_id the start node id
 * @param stop_id the stop node id
 * @param bi_vg the bidirected variation graph
 * @return std::vector<subpaths_t> a vector of paths
 */
subpaths_t get_paths(id_t start_id, id_t stop_id, const bidirected::VariationGraph& bi_vg) {
  std::cout << "[povu::genomics::get_paths]\n";

  // a pair of side and vertex that has been visited
  std::set<std::pair<bidirected::VertexEnd, id_t>> visited;

  // init a visit queue
  bidirected::Vertex const& v_start = bi_vg.get_vertex(start_id);

  const std::set<std::size_t>& edge_from =
	v_start.get_edges_l().size() == 1 ? v_start.get_edges_l() : v_start.get_edges_r();
  const std::set<std::size_t>& edges_to =
	v_start.get_edges_l().size() > 1 ? v_start.get_edges_l() : v_start.get_edges_r();

  std::deque<std::pair<bidirected::VertexEnd, id_t>> q;
  for (auto e_idx: edges_to) {
	const bidirected::Edge& e = bi_vg.get_edge(e_idx);
	e.get_v1_idx() == start_id ?
	  q.push_back(std::make_pair(e.get_v2_end(), e.get_v2_idx())) :
	  q.push_back(std::make_pair(e.get_v1_end(), e.get_v1_idx()));
  }

  std::set<std::pair<bidirected::VertexEnd, id_t>> seen;

  auto unique_push_back_q =
	[&q, &seen, stop_id](std::pair<bidirected::VertexEnd, id_t> e) {
	  if (seen.count(e) || e.second > stop_id) { return; }

	  std::cout << "\t" << "q push back: " << e.first << " " << e.second << "\n";

	  q.push_back(e);
	  seen.insert(e);
	};

  // nodes that are in the path leading to the key vertex
  std::map<side_n_id_t, std::set<side_n_id_t>> in_path;
  for (std::size_t i{start_id}; i < stop_id+1; ++i) {
	in_path[std::make_pair(bidirected::VertexEnd::r, i)] = {};
	in_path[std::make_pair(bidirected::VertexEnd::l, i)] = {};
  }

  bidirected::VertexEnd side =
	v_start.get_edges_l().size() == 1? bidirected::VertexEnd::l : bidirected::VertexEnd::r ;
  in_path[std::make_pair(side, start_id)].insert(std::make_pair(side, start_id));

  // the key is the vertex id and the value is a vector of (unique) paths
  // leading up to the key vertex
  std::map<side_n_id_t, subpaths_t> paths_map;
  // start by adding the start vertex as the only path to itself
  paths_map[std::make_pair(side, start_id)] =
	std::vector<std::vector<side_n_id_t>>{{std::make_pair(side, start_id)}};

  subpaths_t curr_paths;

  auto extend_paths = [&paths_map](side_n_id_t adj,side_n_id_t alt_adj) {
	paths_map[adj] = paths_map[alt_adj];
	for (subpath_t& path_ : paths_map[adj]) {
	  path_.push_back(adj);
	}
  };

  while (!q.empty()) {

	std::pair<bidirected::VertexEnd, id_t> side_id_pair = q.front();
	bidirected::VertexEnd side = side_id_pair.first;
	id_t v_id = side_id_pair.second;
	std::cout << "q.front: (" << side << ", " << v_id << ")\n";

	if (visited.count(side_id_pair)) {
	  q.pop_front();
	  std::cout << "skipping (already visited)\n";
	  continue;
	}

	bidirected::Vertex const& v = bi_vg.get_vertex(v_id);

	bool go_back{false};

	{
	  std::vector<side_n_id_t> to_adjacents;
	  // if l look towards edges is connected to the right, otherwise look towards edges connected to the left
	  // TODO: both sides are r
	  const std::set<std::size_t>& to_adj_edge_idxs = side == bidirected::VertexEnd::l ? v.get_edges_r() : v.get_edges_r();
	  for (auto e_idx: to_adj_edge_idxs) {
		const bidirected::Edge& e = bi_vg.get_edge(e_idx);

		e.get_v1_idx() == v_id ?
		  to_adjacents.push_back(std::make_pair(e.get_v2_end(), e.get_v2_idx())) :
		  to_adjacents.push_back(std::make_pair(e.get_v1_end(), e.get_v1_idx()));
	  }

	  // populate the queue for the next iteration
	  for (side_n_id_t adj: to_adjacents) {
					  unique_push_back_q(adj);

	  }
	}


	const std::set<std::size_t>& frm_adj_edges_idxs = side == bidirected::VertexEnd::l ? v.get_edges_l() : v.get_edges_r();

	std::vector<side_n_id_t> frm_adjacents;
	for (auto e_idx: frm_adj_edges_idxs) {
	  const bidirected::Edge& e = bi_vg.get_edge(e_idx);

	  e.get_v1_idx() == v_id ?
		frm_adjacents.push_back(std::make_pair(e.get_v2_end(), e.get_v2_idx())) :
		frm_adjacents.push_back(std::make_pair(e.get_v1_end(), e.get_v1_idx()));
	}

		// print adjacents
	std::cout << "Frm Adjacents:\n";
	for (auto adj: frm_adjacents) {
	  std::cout << "\t(" << adj.first << ", " << adj.second << ")\n";
	}


	//std::cout << "\tadjacents:\n";
	for (side_n_id_t adj: frm_adjacents) {
	  // std::format("\t\tfrom: {}\n", adj.second);
	  std::cout << "\t\tadj: " << adj.first << ", " << adj.second << "\n";
	  bidirected::VertexEnd adj_side = adj.first;
	  bidirected::VertexEnd alt_adj_side = adj_side == bidirected::VertexEnd::l ? bidirected::VertexEnd::r : bidirected::VertexEnd::l;

	  if (adj.second > stop_id) { continue; }

	  //std::cout << "\t\t" << adj_side << ", " << alt_adj_side << ", " << adj.second << "\n";
	  side_n_id_t alt_adj = std::make_pair(alt_adj_side, adj.second);

	  // print the paths_map keys
	  //std::cout << "\t\tkeys:\n";
	  //for (auto& [k, v]: paths_map) {
	  //std::cout << "\t\t\t" << k.first << ", " << k.second << "\n";
		//}


	  // check whether adj is in paths_map
	  if (paths_map.find(adj) == paths_map.end() && paths_map.find(alt_adj) == paths_map.end() )
	  {
		// push from to queue and break
		std::cout << "\t\tpush front: " << adj.first << ", " << adj.second << "\n";
		q.push_front(adj);
		go_back = true;
		break;
	  }

	  // TODO: move this higher
	  // there is a cycle
	  if (in_path[adj].count(side_id_pair)) {
		std::cout << "\tskipping: cycle detected\n";
		continue;
	  }
	  std::cout << "\t\tnot a cycle\n";


	  if (paths_map.find(adj) == paths_map.end() ) {
extend_paths(adj, alt_adj);
	  }

	  subpaths_t& in_paths = paths_map.at(adj);
	  for (subpath_t& p: in_paths) {
		subpath_t p_ = p;
		p_.push_back(side_id_pair);
		curr_paths.push_back(p_);

		// save that all these vertices are on at least one the paths to v_id
		// populate in_path
				  std::cout << "\t\tinserting:\n";
		for (side_n_id_t x: p_) {
		  std::cout << "\t\t\t" << x.first << " " << x.second << "\n";
		  in_path[side_id_pair].insert(x);
		}

	  }
	}

	if (!go_back) {
	  paths_map[side_id_pair] = curr_paths;
	  visited.insert(side_id_pair);
	}

	curr_paths.clear();
  }

  side = compute_single_edge_side(bi_vg, stop_id);  
  side_n_id_t exit = std::make_pair(side, stop_id);

  // if no path has reached the exit, that is because
  // the opposite side
  if (paths_map[exit].empty()) {
    extend_paths(exit, std::make_pair(complement_side(side), stop_id));
  }

  //std::cout << "returning: " << side << " " << stop_id << "\n";
  subpaths_t simplified_paths = simplify_paths(paths_map[exit]);

  return simplified_paths;
}

/**
 * @brief Given a bidirected::VariationGraph and a subpath, return the string sequence
 */
std::string path_to_seq( const bidirected::VariationGraph& bd_vg, const subpath_t& p) {
  std::string seq = "";
  for (auto [side, id]: p) {
	bidirected::Vertex const& v = bd_vg.get_vertex(id);
	seq += side == bidirected::VertexEnd::l ? v.get_label() : v.get_rc_label();
  }
  
  return seq;
}


// TODO: handle undefined
/**
 * @brief
 */

std::map<std::size_t, std::vector<vcf::vcf_record>>
gen_vcf_records(const bidirected::VariationGraph& bd_vg,
					 const core::config& app_config,
					 const std::vector<std::vector<std::set<std::size_t>>>& haplotypes,
					 const std::vector<subpaths_t>& all_paths) {

  std::vector<std::string> reference_paths = app_config.get_reference_paths();
  std::map<std::size_t, std::vector<vcf::vcf_record>> vcf_records;


  std::map<std::string, std::size_t> paths_map; // path to id
  std::map<std::size_t, std::string> paths_map_rev; //  id to path name
  for (auto p : bd_vg.get_paths()) {
	paths_map[p.name] = p.id;
	paths_map_rev[p.id] = p.name;
  }


  for (auto ref_name : reference_paths) {
	std::size_t ref_id = paths_map[ref_name];
	//std::cout << "ref_id: " << ref_id << " ref name: " << ref_name << "\n";
	vcf_records[ref_id] = std::vector<vcf::vcf_record>();
  }



  if (true) {
	// add the undefined reference for flubbles that don't have the a reference
	paths_map[UNDEFINED_PATH_LABEL] = UNDEFINED_PATH_ID;
		paths_map_rev[UNDEFINED_PATH_ID] = UNDEFINED_PATH_LABEL;
	reference_paths.push_back(UNDEFINED_PATH_LABEL);
  }



  std::size_t num_flubbles = haplotypes.size();
  std::size_t num_paths{};
  // for each flubble
  for (std::size_t fl_idx{}; fl_idx < num_flubbles; ++fl_idx) {
	num_paths = all_paths[fl_idx].size();

	std::cout << "flubble: " << fl_idx << " num paths " << num_paths << "\n";

	// does this flubble contain a reference path?
	bool has_ref{false};

	// while being strand aware
	// spell all the sequences in the flubble as a vector of strings
	// where each string is a path
	std::vector<std::string> path_seqs;
	for (std::size_t path_idx{}; path_idx < num_paths; ++path_idx) {
	  const subpath_t& p = all_paths[fl_idx][path_idx];
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

	  // every path should be supported by at least one haplotype
	  if (path_haplotypes.empty()) {
		std::cerr << "ERROR: path " << path_idx << " in flubble " << fl_idx << " has no haplotypes\n";
	  }

	  // because we want a VCF for each reference path then
	  // for each reference path we check if it is in the set of haplotypes
	  // for that path
	  // if it is then remove it from the path seqs and the rest are variants
	  // we then populate the vcf records map which contains
	  // the vcf records for each reference path
	  for (const std::string& ref_name : reference_paths) {
		std::size_t ref_id = paths_map[ref_name];
		if (path_haplotypes.count(ref_id)) {

		  has_ref = true;
		  
		  // the path contains the reference haplotype ref_id
		  vcf::vcf_record vcf_rec;

		  // get the position of that path in the reference whose id is ref_id
		  id_t v_id = all_paths[fl_idx][path_idx].front().second;
		  bidirected::Vertex const& v = bd_vg.get_vertex(v_id);
		  // TODO: can this loop be avoided?
		  for (const bidirected::PathInfo& vertex_path : v.get_paths()) {
			if (vertex_path.path_id == ref_id) {
			  vcf_rec.pos = vertex_path.step_index;
			  break;
			}

			vcf_rec.ref = path_seqs[path_idx];
			vcf_rec.alt = immutable_erase(path_seqs, path_idx);

			vcf_records[ref_id].push_back(vcf_rec);
		  }
		}
	  }	  
	}
	

	if (!has_ref) {
	  vcf::vcf_record vcf_rec;
	  vcf_rec.pos = UNDEFINED_PATH_POS;
	  vcf_rec.ref = "";
	  vcf_rec.alt = path_seqs;


	  vcf_records[UNDEFINED_PATH_ID].push_back(vcf_rec);
	}
  }
  return vcf_records;
}

/**
 * Associate each path with a set of haplotypes
 *
 */
std::vector<std::vector<std::set<std::size_t>>>
find_path_haplotypes(const std::vector<subpaths_t>& all_paths, const bidirected::VariationGraph& bd_vg) {
  std::vector<std::vector<std::set<std::size_t>>> haplotypes_per_path;

  for (const subpaths_t& c_flubble_paths: all_paths) {
	std::vector<std::set<std::size_t>> h;

	for (const subpath_t& c_flubble_path: c_flubble_paths) {

	  // for each path in the flubble
	  // find the haplotype
	  // add the haplotype to the path
	  std::set<std::size_t> curr_haplotypes;
	  std::set<std::size_t> intersect;
	  std::set<std::size_t> temp;
	  std::vector<std::vector<std::size_t>> haplotypes_per_path;
	  for (const side_n_id_t& side_id_pair: c_flubble_path) {
		// find the haplotype
		const std::vector<bidirected::PathInfo>& path_info =
		  bd_vg.get_vertex(side_id_pair.second).get_paths();

		for (auto p: path_info) {
		  curr_haplotypes.insert(p.path_id);
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



void call_variants(const tree::Tree& pvst_, const bidirected::VariationGraph& bd_vg, const core::config& app_config) {
  if (app_config.verbosity() > 3) {
	std::cout << "[povu::genomics::call_variants]" << std::endl;
  }

  //output_format of = output_format::VCF;
  std::vector<std::pair<std::size_t, std::size_t>> canonical_flubbles =
	extract_canonical_flubbles(pvst_);

  // walk paths in the digraph
  // while looping over canonical_flubbles

  std::vector<subpaths_t> all_paths;

  // extract flubble paths
  for (auto c_flubble: canonical_flubbles) {
	if (app_config.verbosity() > 4) {
	  std::cout << std::format("\tcanonical flubble: ({}, {})\n", c_flubble.first, c_flubble.second);
	}

	// std::size_t offset = 1;
	// TODO: remove these assignments
	std::size_t start = c_flubble.first;
	std::size_t stop = c_flubble.second;

	// limited DFS
	subpaths_t paths = get_paths(start, stop, bd_vg);

	{

	  std::cout << "\tNumber of paths: " << paths.size() << std::endl;

	  // print paths
	  if (app_config.verbosity() > 4) {
		std::cout << "\t";
		for (auto p: paths) {
		  for (side_n_id_t x: p) {
			std::cout << x.first << " " << x.second << ",  ";
		  }
		  std::cout << "\n";
		}
		std::cout << std::endl;
	  }
	}

	all_paths.push_back(paths);
  }

  std::vector<std::vector<std::set<std::size_t>>> haplotypes_per_path =
	find_path_haplotypes(all_paths, bd_vg);

  std::map<std::size_t, std::vector<vcf::vcf_record>> vcf_records =
	gen_vcf_records(bd_vg, app_config, haplotypes_per_path, all_paths);

  vcf::write_vcfs(vcf_records, bd_vg);


  {
  /*
  // TODO: should this be in digraph?

  // include the undefined reference

   std::size_t ref_id{};

   for (std::string ref_path : app_config.reference_paths) {

	 std::cout << "\n\n";
	 //std::cout << "ref_path: " << ref_path << std::endl;

	ref_id =  paths_map[ref_path];
	// TODO: move to print VCF
	vcf::print_vcf_header(ref_path);

	for (std::size_t i{0}; i < all_paths.size(); ++i) {
	  std::vector<std::vector<std::size_t>> paths = all_paths[i];
	  bool print_res  =
		vcf::print_vcf(paths, dg, ref_path, ref_id);

	  if (print_res) {
		// delete the path from all_paths
		all_paths.erase(all_paths.begin()+i);
		--i;
	  }
	}

	//for (auto paths: all_paths) {
	//  vcf::print_vcf(paths, dg, ref_path, ref_id);
	//}
  }
  */
	}
}

} // namespace genomics
