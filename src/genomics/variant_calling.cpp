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

#include "../pvst/pvst.hpp"
#include "../graph/tree.hpp"
#include "../graph/digraph.hpp"
#include "../core/constants.hpp"
#include "../core/core.hpp"
#include "./genomics.hpp"

namespace genomics {

// TODO: this can be done while constructing the PVST
std::vector<std::pair<std::size_t, std::size_t>>
extract_canonical_flubbles(tree::Tree pvst_) {

  tree::Vertex const& root = pvst_.get_root();

  // subtree set

  std::vector<std::pair<std::size_t, std::size_t>> canonical_flubbles;
  //std::vector<std::size_t> canonical_subtrees;
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
	  //canonical_subtrees.push_back(current_vertex);
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

void call_variants(
  tree::Tree pvst_,
  digraph::DiGraph dg,
  core::config app_config) {

  std::cout << "[genomics::call_variants]" << std::endl;
    
  //output_format of = output_format::VCF;
  std::vector<std::pair<std::size_t, std::size_t>> canonical_flubbles =
	extract_canonical_flubbles(pvst_);

  // walk paths in the digraph
  // while looping over canonical_flubles
  
  std::vector<std::vector<std::vector<std::size_t>>> all_paths;

  // extract flubble paths
  for (auto c_flubble: canonical_flubbles) {

	//std::cout << "canonical flubble: " << c_flubble.first << " " << c_flubble.second << std::endl;
	
	std::size_t offset = 1;
	std::size_t start = c_flubble.first-offset;
	std::size_t stop = c_flubble.second - offset;

	//limited_dfs(start, stop, dg);
	std::vector<std::vector<std::size_t>> paths = get_paths(start, stop, dg);
	all_paths.push_back(paths);

  }

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
}

} // namespace genomics

