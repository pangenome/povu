#include <cstddef>
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


#include "../pvst/pvst.hpp"
#include "../graph/tree.hpp"
#include "../graph/digraph.hpp"
#include "../core/constants.hpp"
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

	std::cout << "current vertex: " << current_vertex << std::endl;

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
  std::cout << "Canonical subtrees:" << std::endl;
  for (auto subtree: canonical_flubbles) {
	std::cout << subtree.first << " " << subtree.second << std::endl;
  }

  return canonical_flubbles;
  // use canonical subtrees to make a vector of pairs of flubble starts and stops
  //std::vector<std::pair<std::size_t, std::size_t>> canonical_flubbles;
}


std::vector<std::vector<std::size_t>>
get_paths(std::size_t start, std::size_t stop, digraph::DiGraph dg) {

  // typedef std::vector<std::vector<std::size_t>> x;
  std::vector<std::vector<std::size_t>> paths;

  // the paths leading up to key vertex
  std::map<std::size_t, std::vector<std::vector<std::size_t>>> paths_map;

  // nodes that are in the path leading to the key vertex
  std::map<std::size_t, std::set<std::size_t>> in_path;
  for (std::size_t i{start}; i < stop+1; ++i) {
	in_path[i] = {};
  }
  

  std::vector<std::size_t> path {9};
  //path.push_back(start);
  paths.push_back(path);

  paths_map[start] = paths;

  in_path[start].insert(start);
  
  for (std::size_t i{start+1}; i < stop+1; ++i) {
	// get the vertex at index i

	//std::cout << "\t" << "i: " << i << "\n";

	digraph::Vertex const& v = dg.get_vertex(i);
	  std::vector<std::vector<std::size_t>>	curr_paths;
	  
	for (auto adj: v.in()) {
	  //std::cout << "\t\t" << "frm:" << adj.from() << "\n";
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


  
void call_variants(tree::Tree pvst_, digraph::DiGraph dg) {
  //output_format of = output_format::VCF;
  std::vector<std::pair<std::size_t, std::size_t>> canonical_flubbles =
	extract_canonical_flubbles(pvst_);

  // walk paths in the digraph
  // while looping over canonical_flubles

  vcf::print_vcf_header();
  
  for (auto c_flubble: canonical_flubbles) {
	std::size_t offset = 1;
	std::size_t start = c_flubble.first-offset;
	std::size_t stop = c_flubble.second - offset;

	//limited_dfs(start, stop, dg);
	std::vector<std::vector<std::size_t>> paths = get_paths(start, stop, dg);


	vcf::print_vcf_variants(paths, dg);
	//v.out();
  }

}


} // namespace genomics

