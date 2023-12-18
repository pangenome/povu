#include <algorithm>
#include <bit>
#include <cstddef>
#include <iostream>
#include <map>
#include <string>
#include <sys/types.h>
#include <vector>
#include <list>

#include "./pvst.hpp"
#include "../graph/tree.hpp"
#include "../core/core.hpp"

namespace pvst {
const std::size_t INVALID_ID = core::constants::UNDEFINED_SIZE_T;

/**
 * @brief compute_pst
 *
 *  The input is a vector of cycle equivalence classes.
 *  The output is a PST
 *
 * @param v[in] vector of cycle equivalence classes
 * @return tree::Tree the PST
*/
tree::Tree compute_pst(std::vector<std::size_t> classes) {
  // TODO this tree is lager than should be
  tree::Tree t = tree::Tree();

  // get the equivalence class of the ith element
  auto get_class = [&classes](std::size_t i) -> std::size_t {
	return classes[i];
  };

  // a map of equivalence classes to a vector of their positions in the edge stack
  std::map<size_t, std::vector<size_t>> pos_map;
  for (std::size_t i{}; i < classes.size(); ++i) {
	pos_map[get_class(i)].push_back(i);
  }

  std::size_t counter{}; // the number of nodes in the tree
  std::size_t curr_idx{}; // the index of the current node in the tree
  std::size_t curr_class{}; // the current equivalence class
  std::size_t prev_class{}; // the previous equivalence class

  // a vector of bools with a t or f on whether a class is nesting or not
  // values are true when we are in a SESE region
  std::vector<bool> nesting(classes.size(), false);

  // the value of pos map at the current index
  std::vector<std::size_t> temp{};
  // get the current position in the edge stack of the current equivalence class
  // e.g this is the nth - 1 time we are seeing the current equivalence class
  // i is an index in v
  auto find_pos = [&](std::size_t i) -> std::size_t {
	for (std::size_t j{}; j < temp.size(); ++j) {
	  if (i == temp[j]) { return j; }
	}
	// TODO: what happens in this case? throw an exception?
	return core::constants::UNDEFINED_SIZE_T;
  };

  // this is the nth - 1 time we are seeing the current equivalence class
  std::size_t pos;

  // the nth - 1 time we saw the previous equivalence class
  std::size_t prev_pos;

  std::size_t parent_idx{}; // TODO: replace with p

  // a lambda that takes a parent index and updates the tree
  auto update_tree = [&](std::size_t parent) -> void {
	t.add_vertex(parent, counter, curr_class);
	curr_idx = counter;
	++counter;
  };

  for (std::size_t i{}; i < classes.size(); ++i) {
	curr_class = get_class(i);

	temp = pos_map[curr_class];

	// if the class only has one element, then it is not a SESE region
	if (temp.size() == 1) { continue; }

	pos = find_pos(i);

	if (t.empty()) {
	  t.add_vertex(core::constants::UNDEFINED_SIZE_T, counter, curr_class);
	  ++counter;
	}
	else if (nesting[curr_class]) {
	  // if we are in a SESE region
	  // we are at the end of the SESE region
	  // When a region is exited, the current region is set to be the exited region’s parent.
	  // curr_idx has class curr_class
	  curr_idx = t.get_parent(curr_idx);
	  nesting[curr_class] = false;
	}
	else if (prev_class == curr_class) {
	  // if the same class occurs in tandem
	  // then we are in a SESE region

	  // this is not the last time we are seeing the current class
	  // if the current class occurs again later
	  // i.e pos is not the last value in temp
	  //
	  // or
	  //
	  // we expect the last node to have the same class as the current class if not add it because we
	  // have two nodes with the same class after:
	  //  - a nesting region 2 ... 2 2
	  if (pos < temp.size() - 1
		|| t.get_vertex(counter -1).get_class() != curr_class) {
		std::size_t p = t.get_parent(curr_idx);
		update_tree(p);
	  }
	}
	else if ((prev_class != curr_class))  {
	  // we may entering a new SESE region
	  // whereby prev_class is the parent of curr_class

	  // this is not the last time we are seeing the previous class
	  if (prev_pos < pos_map[prev_class].size() - 1 ) {
		nesting[prev_class] = true;
		update_tree(curr_idx);
	  }
	  else if (pos < temp.size() - 1) {
		std::size_t p = t.get_parent(curr_idx);
		update_tree(p);
	  }
	  else {
		// that is the last time we are seeing the previous class
		// so we are exiting the SESE region
		curr_idx = t.get_parent(curr_idx);
	  }
	}
	// assign previous class from the current class
	prev_class = curr_class;
	prev_pos = pos;
  }

  return t;
}

/**
 * @brief compute_pvst
 *
 * A PVST a tree in which each subtree is a flubble
 * The input is a vector of tuples. Each tuple contains the following information:
 *   - the first element is the eq class
 *   - the second element is the one node
 *   - the third element is the other edge
 *
 *
 * @param v[in] vector of tuples
 * @return tree::Tree the PST
*/

tree::Tree compute_pvst(std::vector<std::tuple<std::size_t, std::size_t, std::size_t>> const& v) {
  // TODO this tree is lager than should be
  tree::Tree t = tree::Tree();

  // get the equivalence class of the ith element
  auto get_class = [&v](std::size_t i) -> std::size_t {
	return std::get<0>(v[i]);
  };

  // get_src
  auto get_src = [&v](std::size_t i) -> std::size_t {
	return std::get<1>(v[i]);
  };

  // get tgt
  auto get_tgt = [&v](std::size_t i) -> std::size_t {
	return std::get<2>(v[i]);
  };

  // a map of equivalence classes to a vector of their positions in the edge stack
  std::map<size_t, std::vector<size_t>> pos_map;
  for (std::size_t i{}; i < v.size(); ++i) {
	pos_map[get_class(i)].push_back(i);
  }

  std::size_t counter{}; // the number of nodes in the tree
  // curr idx is the last added node in the tree or if we have exited an
  // flubble it is the parent of the last added node
  std::size_t curr_idx{}; // the index of the current node in the tree
  std::size_t curr_class{}; // the current equivalence class
  std::size_t prev_class{}; // the previous equivalence class

  // a vector of bools with a t or f on whether a class is nesting or not
  // values are true when we are in a SESE region
  std::vector<bool> nesting(v.size(), false);

  // the value of pos map at the current index
  std::vector<std::size_t> temp{};
  // get the current position in the edge stack of the current equivalence class
  // e.g this is the nth - 1 time we are seeing the current equivalence class
  // i is an index in v
  auto find_pos = [&](std::size_t i) -> std::size_t {
	for (std::size_t j{}; j < temp.size(); ++j) {
	  if (i == temp[j]) { return j; }
	}
	// TODO: what happens in this case? throw an exception?
	return core::constants::UNDEFINED_SIZE_T;
  };

  // this is the nth - 1 time we are seeing the current equivalence class
  std::size_t pos;

  // the nth - 1 time we saw the previous equivalence class
  std::size_t prev_pos;

  // the parent index of the current node, not always pointing to the parent of
  // the current node
  std::size_t parent_idx{};

  // source and target of the current edge as strings
  std::string src;
  std::string tgt;

  std::set<std::size_t> seen;

  // a lambda that takes a parent index and updates the tree
  auto update_tree = [&](std::size_t parent, std::string& meta, bool update_curr_idx = true) -> void {
	// convert meta to a size_t and if check if it has been seen
	// if it has been seen then we do not add it to the tree
	// if it has not been seen then we add it to the tree
	std::size_t meta_size_t = std::stoull(meta);
	if (seen.find(meta_size_t) != seen.end()) { return; }
	seen.insert(meta_size_t);

	t.add_vertex(parent, counter, curr_class, meta);


	// swap meta values between parent and counter nodes
	if (std::stoull(t.get_vertex(parent).get_meta()) > std::stoull(t.get_vertex(counter).get_meta())) {
	  std::string temp = t.get_vertex(parent).get_meta();
	  t.get_vertex_mut(parent).set_meta(t.get_vertex(counter).get_meta());
	  t.get_vertex_mut(counter).set_meta(std::move(temp));
	}


	if (update_curr_idx) { curr_idx = counter; }
	++counter;
  };

  for (std::size_t i{}; i < v.size(); ++i) {
	curr_class = get_class(i);
	src = std::to_string(get_src(i));
	tgt = std::to_string(get_tgt(i));

	temp = pos_map[curr_class];

	pos = find_pos(i);

	if (t.empty()) {
	  t.add_vertex(core::constants::UNDEFINED_SIZE_T, counter, curr_class,src);
	  ++counter;

	  std::size_t meta_size_t = std::stoull(src);
	  seen.insert(meta_size_t);
	}
	else if (nesting[curr_class]) {
	  // if we are in a SESE region
	  // we are at the end of the SESE region
	  // When a region is exited, the current region is set to be the exited region’s parent.
	  // curr_idx has class curr_class
	  curr_idx = t.get_parent(curr_idx);
	  update_tree(curr_idx, src, false);
	  nesting[curr_class] = false;
	}
	else if (prev_class == curr_class) {
	  //
	  parent_idx = t.get_parent(curr_idx);
	  update_tree(parent_idx, tgt);
	  update_tree(parent_idx, src);

	}
	else if (prev_class != curr_class) {
	  // we may entering a new SESE region
	  // whereby prev_class is the parent of curr_class

	  // add the node of the current class to the tree as a child of the previous class

	  // this is not the last time we are seeing the prev class
	  if (prev_pos < pos_map[prev_class].size() - 1 ) {
		nesting[prev_class] = true;
		update_tree(curr_idx, src, false);
		update_tree(curr_idx, tgt);
	  }
	  // this is not the last time we are seeing the current class
	  //else if (pos < temp.size() - 1) { }
	  // the previous class is not the parent of the current class
	  // rather they are just adjacent
	  // e.g 9 9 10 10 or 8 9
	  else {
		parent_idx = t.get_parent(curr_idx);
		update_tree(parent_idx, tgt);

		// that is the last time we are seeing the previous class
		// so we are exiting the SESE region

		// should be handled by the nesting case
		//curr_idx = t.get_parent(curr_idx);
	  }
	}

	// assign previous class from the current class
	prev_class = curr_class;
	prev_pos = pos;
  }

  return t;
}

tree::Tree compute_pvst(std::vector<std::pair<std::size_t, std::size_t>> v, const core::config& app_config)   {
  std::string fn_name = "[povu::pvst::compute_pvst]";

  std::vector<core::eq_n_id_t> v_;
  std::transform(v.begin(), v.end(),
				v_.begin(), [](auto& p) -> core::eq_n_id_t { return {p.second, p.first}; });
  return compute_pvst(v_, app_config);
}



/**
 * @brief compute the pvst of a given vector of eq_n_id_t values
 *
 * @param v a vector of pairs
 *          where each pair is (index, eq class)
 * @return a tree
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

	if (!t.empty()) {
	  std::cout << "root: " << t.get_root().get_id() << std::endl;
	}
	
	if (true) {
	  std::cerr << "i: " << i

				<< " counter " << counter
						<< " p_id: " << p_id
		//<< " be idx: " << be_idx
		//<< " bd idx: " << bd_idx
				<< " current class: " << current_class

		//<< " bd idx str: " << bd_idx_str
				<< std::endl;
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
	  // in this case the eq class nests something
	  std::list<std::size_t> positions;
	  std::copy(ps.begin(), ps.end(), std::back_inserter(positions));

	  for (auto it = positions.begin(); it != positions.end(); ++it) {

		// detect a bubble chain
		// handle bubble chains
		if (*it == i && std::next(it) != positions.end() && it != positions.begin()
			//&& (*it  - *std::prev(it ) > 1)
			&& ((*std::next(it) - *it ) > 1) ) {
		  // add a dummy parent?

		  tree::Vertex& prt_v = t.get_vertex_mut(p_id);
		  
		  if (!prt_v.is_dummy()) {
			
			if (is_root(p_id)) {
			  std::string m = prt_v.get_meta();
			  ++counter;
			  t.add_vertex(INVALID_ID, counter, prt_v.get_class(), m, true);
			  t.set_root(counter);
			  std::cout << "ctr" << counter << std::endl;
			t.get_vertex_mut(p_id).set_parent(counter);
			  prt_v.set_parent(counter);
			  std::cout << prt_v.get_parent() << " p " << t.get_vertex(p_id).get_parent() << std::endl;
			  
			  t.get_vertex_mut(counter).add_child(p_id);
			  t.add_vertex(p_id, counter, current_class, bd_idx_str, true);

			  p_class = prt_v.get_class();
			}
			else {
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

	
	
	if ((!nesting[current_class] && ps.size() > 1) ||
		(!nesting[current_class] && i+1 < v.size() && current_class != v[i+1].eq_class)) {
	  std::cout << "down\n";
	  nesting[current_class] = true;
	  p_class = current_class;
	  p_id = counter;
	}

	
	if (ps.back() == i || (!t.get_vertex(counter).is_dummy() && current_class == p_class && i+1 < v.size() && current_class != v[i+1].eq_class)) {
	  	  	  std::cout << "flip\n";
	  nesting[current_class] = false;
	}

	if (!nesting[p_class] && !is_root(p_id)) {
	  	  std::cout << "up\n";
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
} // namespace pvst
