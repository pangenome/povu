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


// TODO: this can be done while constructing the PVST
std::vector<std::pair<std::size_t, std::size_t>>
extract_canonical_flubbles(tree::Tree pvst_) {

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

	if (children.empty()) { continue; }

	is_canonical_subtree = true;

	max_child = core::constants::SIZE_T_MIN;

	for (auto child: children) {
	  std::set<std::size_t> const& c =  pvst_.get_children(child);

	  if (child > max_child) { max_child = child; }

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

std::vector<std::vector<std::size_t>>
get_paths_old(std::size_t start, std::size_t stop, digraph::DiGraph dg) {
  std::cout << "[genomics::get_paths]\n";

  //std::cout << "\t" << "start: " << start << " stop " << stop << std::endl;

  // typedef std::vector<std::vector<std::size_t>> x;
  std::vector<std::vector<std::size_t>> paths;

  // the paths leading up to key vertex
  std::map<std::size_t, std::vector<std::vector<std::size_t>>> paths_map;

  // nodes that are in the path leading to the key vertex
  std::map<std::size_t, std::set<std::size_t>> in_path;
  for (std::size_t i{start}; i < stop+1; ++i) {
	in_path[i] = {};
  }


  std::vector<std::size_t> path {start};
  //path.push_back(start);
  paths.push_back(path);

  paths_map[start] = paths;

  in_path[start].insert(start);

  for (std::size_t i{start+1}; i < stop+1; ++i) {
	// get the vertex at index i

	std::cout << "\t" << "i: " << i << "\n";

	digraph::Vertex const& v = dg.get_vertex(i);
	  std::vector<std::vector<std::size_t>>	curr_paths;

	for (auto adj: v.in()) {
	  std::cout << "\t\t" << "frm:" << adj.from() << "\n";
	  std::vector<std::vector<std::size_t>>& in_paths = paths_map.at(adj.from());

	  // there is a cycle
	  if (in_path[adj.from()].count(i) ) { continue; }

	  for (std::vector<std::size_t>& p: in_paths) {
		std::vector<std::size_t> p_ = p;
		//curr_paths.push_back
		p_.push_back(i);
		curr_paths.push_back(p_);

		// populate in_path
		for (auto v: p_) {
		  in_path[i].insert(v);
		}

	  }
	  /*
	  std::cout << "\t\t" << "curr_paths: " << "\n" << "\t\t\t";
	  for (auto p: curr_paths) {
		for (auto v: p) {
		  std::cout << v << " ";
		}
		std::cout << "\n";
		}

	  */
	}
	paths_map[i] = curr_paths;

  }

  // print the content of paths_map at stop

  /*
  std::cout << "Final:\n";
  for (auto pp: paths_map[stop]) {
	for (auto p: pp) {
	  std::cout << p << " ";
	}
	std::cout << "\n";
  }
  */


  return paths_map[stop];
}

std::vector<std::vector<std::size_t>>
get_paths(std::size_t start, std::size_t stop, digraph::DiGraph dg) {
  std::cout << "[povu::genomics::get_paths]\n";

  //std::cout << "\t" << "start: " << start << " stop " << stop << std::endl;

  // typedef std::vector<std::vector<std::size_t>> x;
  std::vector<std::vector<std::size_t>> paths;

  // the paths leading up to key vertex
  std::map<std::size_t, std::vector<std::vector<std::size_t>>> paths_map;

  // nodes that are in the path leading to the key vertex
  std::map<std::size_t, std::set<std::size_t>> in_path;
  for (std::size_t i{start}; i < stop+1; ++i) { in_path[i] = {}; }


  std::vector<std::size_t> path {start};
  //path.push_back(start);
  paths.push_back(path);

  paths_map[start] = paths;

  in_path[start].insert(start);

  // push every vertex into a queue
  // a queue of size_t
  std::deque<std::size_t> q;
  for (std::size_t i{start+1}; i < stop+1; ++i) {
	q.push_back(i);
  }

  std::set<std::size_t> visited;

  while (!q.empty()) {
	std::size_t i = q.front();

	if (visited.count(i)) {
	  q.pop_front();
	  continue;
	}

	// get the vertex at index i
	std::cout << "\t" << "i: " << i << "\n";

	digraph::Vertex const& v = dg.get_vertex(i);

	std::vector<std::vector<std::size_t>> curr_paths;

	bool go_back{false};

	for (auto adj: v.in()) {
	  std::cout << "\t\t" << "frm:" << adj.from() << "\n";

	  // check whether adj.from() is in paths_map
	  if (paths_map.find(adj.from()) == paths_map.end() ) {
		// push from to queue and break
		q.push_front(adj.from());
		go_back = true;
		break;
	  }

	  std::vector<std::vector<std::size_t>>& in_paths = paths_map.at(adj.from());

	  // there is a cycle
	  if (in_path[adj.from()].count(i) ) { continue; }

	  for (std::vector<std::size_t>& p: in_paths) {
		std::vector<std::size_t> p_ = p;
		//curr_paths.push_back
		p_.push_back(i);
		curr_paths.push_back(p_);

		// populate in_path
		for (auto v: p_) {
		  in_path[i].insert(v);
		}

	  }
	  /*
	  std::cout << "\t\t" << "curr_paths: " << "\n" << "\t\t\t";
	  for (auto p: curr_paths) {
		for (auto v: p) {
		  std::cout << v << " ";
		}
		std::cout << "\n";
		}

	  */
	}

	if (!go_back) {
	  paths_map[i] = curr_paths;
	  visited.insert(i);
	}


  }


  // print the content of paths_map at stop

  /*
  std::cout << "Final:\n";
  for (auto pp: paths_map[stop]) {
	for (auto p: pp) {
	  std::cout << p << " ";
	}
	std::cout << "\n";
  }
  */


  return paths_map[stop];
}

subpaths_t get_paths(id_t start_id, id_t stop_id, bidirected::VariationGraph bi_vg) {
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

  auto complement_side = [](bidirected::VertexEnd side) {
	return side == bidirected::VertexEnd::l ? bidirected::VertexEnd::r : bidirected::VertexEnd::l;
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

  bidirected::Vertex const& v_stop = bi_vg.get_vertex(stop_id);
  side = v_stop.get_edges_l().size() == 1? bidirected::VertexEnd::l : bidirected::VertexEnd::r;

  side_n_id_t exit = std::make_pair(side, stop_id);
  
  if (paths_map[exit].empty()) {
	extend_paths(exit, std::make_pair(complement_side(side), stop_id));
  }
  
  std::cout << "returning: " << side << " " << stop_id << "\n";
  
 return paths_map[exit];
}

void call_variants(tree::Tree pvst_, bidirected::VariationGraph bd_vg, core::config app_config) {
  if (app_config.verbosity() > 3) {
	std::cout << "[povu::genomics::call_variants]" << std::endl;
  }

  //output_format of = output_format::VCF;
  std::vector<std::pair<std::size_t, std::size_t>> canonical_flubbles =
	extract_canonical_flubbles(pvst_);

  // walk paths in the digraph
  // while looping over canonical_flubles

  std::vector<subpaths_t> all_paths;

  // extract flubble paths
  for (auto c_flubble: canonical_flubbles) {
	if (app_config.verbosity() > 4) {
	  std::cout << std::format("\tcanonical flubble: ({}, {})\n",
							   c_flubble.first, c_flubble.second);
	}

	// std::size_t offset = 1;
	// TODO: remove these assignments
	std::size_t start = c_flubble.first;
	std::size_t stop = c_flubble.second;

	// limited DFS
	subpaths_t paths = get_paths(start, stop, bd_vg);

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

	all_paths.push_back(paths);
  }

  {
  /*
  // TODO: should this be in digraph?
  std::map<std::string, std::size_t> paths_map; // path to id

  for (auto p : dg.get_paths()) { paths_map[p.name] = p.id; }

  // include the undefined reference
  if (true) {
	// add the undefined reference for flubbles that don't have the a reference
	paths_map["undefined"] = core::constants::UNDEFINED_SIZE_T;
	app_config.reference_paths.push_back("undefined");
  }

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
