#include <bit>
#include <cstddef>
#include <iostream>
#include <map>
#include <string>
#include <vector>
#include <list>

#include "../graph/tree.hpp"
#include "./pvst.hpp"
#include "../core/core.hpp"



namespace pvst {
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
tree::Tree compute_pvst(
  std::vector<std::tuple<std::size_t , std::size_t, std::size_t>> const& v) {
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

/**
 * @brief compute the pvst of a given vector of pairs
 *
 * @param v a vector of pairs
 *          where each pair is (index, eq class)  
 * @return a tree
 */
tree::Tree compute_pvst(std::vector<std::pair<std::size_t, std::size_t>> v) {
  if (true) {
	std::cout << "[povu::pvst::compute_pvst]\n";
  }
  
  tree::Tree t = tree::Tree();

  // lambda to map back an index from the biedged and dummies vertex idx to the
  // bidirected idx
  auto bidirected_idx = [](std::size_t x) -> std::size_t {
	x = x - 1;
	return x % 2 == 0 ? ((x + 2) / 2) - 1 : ((x + 1) / 2) - 1;
  };

	
  for (auto i: v) {
	std::cout << i.first << " b " << bidirected_idx(i.first) << " " << i.second << std::endl; 
  }

  // a map of equivalence classes to a vector of their positions in the edge stack
  std::map<size_t, std::vector<size_t>> pos_map;
  for (std::size_t i{}; i < v.size(); ++i) {
	pos_map[v[i].second].push_back(i);
  }

  // print pos_map
  std::cout << "pos_map\n";
  for (auto i: pos_map) {
	std::cout << i.first << ": ";
	for (auto j: i.second) { std::cout << j << ", "; }
	std::cout << std::endl;
  }

  std::string bd_idx_str;

  std::size_t
	current_parent_class{},
	current_parent_id{},
	current_class{},
	prev_class{},
	bd_idx{}, // bidirected index
	be_idx{}  // biedged index
	;

  tree::Vertex p_v;
  
  std::size_t counter{}, i{};

  while (i < v.size()) {


	
	be_idx = v[i].first;
    bd_idx = bidirected_idx(be_idx);
	current_class = v[i].second;
	bd_idx_str = std::to_string(bd_idx+1);


	std::cout << "i: " << i << " counter: " << counter
			  << " cpid: " << current_parent_id << " be idx " << be_idx << " bd idx: " << bd_idx
			  << std::endl;
	
	// ----------------------
	// handle the root vertex
	// ----------------------

	//
	{
	  if (i==0 && t.empty()) {
		t.add_vertex(core::constants::UNDEFINED_SIZE_T, counter, current_class, bd_idx_str);
		current_parent_class = current_class;
		current_parent_id = i;
		++counter;
		++i;
		continue;
	  }
	  else if (i==0 && !t.empty()) {
		throw std::runtime_error("adding vertices to a rootless tree");
	  }
	}

	// ------------------
	// add other vertices
	// ------------------
	{
	  t.add_vertex(current_parent_id, counter, current_class, bd_idx_str);
	}	

    // ------------------
	// update/determine the parent
	// ------------------


	std::vector<std::size_t> const& ps = pos_map[current_class];
	
	std::list<std::size_t> positions;
	for (auto p: ps) {
	  positions.push_back(p);
	}
	
	if (positions.size () == 1) {
	  ++counter;
	  ++i;
	  continue; }

	bool foo{false};

	for (auto it = positions.begin(); it != positions.end(); ++it) {
	  if (std::next(it) == positions.end()) { break; }

	  if (*it == i && it != positions.begin() &&
		  (*it  - *std::prev(it ) > 1) &&
		((*std::next(it) - *it ) > 1) ) {
		// add a dummy parent?

		auto prt_v = t.get_vertex(current_parent_id);

		if (!prt_v.is_dummy()) {

		  //bool h{false};
		  std::vector<std::size_t> children_found;
		  for (auto c : prt_v.get_children()) {
			// check whether the string contains _d
		  
			if (!t.get_vertex(c).is_dummy() ) {
			  //h = true;
			  ++counter;

			  std::string dummy_str = std::to_string(prt_v.get_id() + 1);
			  t.add_vertex(current_parent_id, counter,
						   prt_v.get_class(), dummy_str, true);
			  break;
			}
		  }

		  for (auto c : prt_v.get_children()) {
			// check whether the string contains _d
		  
			if (!t.get_vertex(c).is_dummy()) {
			  t.get_vertex_mut(c).set_parent(counter);
			  t.get_vertex_mut(counter).add_child(c);
			  children_found.push_back(c);
			  //prt_v.remove_child(c);
			}
		  }

		  for (auto c : children_found) {
			//std::cout << "\t remove "<< c << "\n";
			t.get_vertex_mut(current_parent_id).remove_child(c);
		  }
		}

		// move up the tree to find a non dummy parent
		while (t.get_vertex(current_parent_id).is_dummy() ) {
		  current_parent_id = t.get_parent(current_parent_id);
		}
		
		std::cout << "\tcpid "<< current_parent_id << " " << counter<< "\n";
		++counter;
		std::string dummy_str = std::to_string( (bidirected_idx(v[i].first) + 1));

		t.add_vertex(current_parent_id, counter, current_class, dummy_str, true);

		current_parent_id = counter;
		++i;
		++counter;
		foo = true;
		// update all the 
	  }
	}

	if (foo) {
	  std::cout << "foo was true" << "\n";
	  continue;
	}
	
	if (positions.back() == i) {
		// this is the last time we are seeing this class
		std::cout << "up " << i << "\n";
		// go up
		std::size_t parent_id = t.get_parent(current_parent_id);
		p_v = t.get_vertex(parent_id);
		
		current_parent_class = p_v.get_class();
		current_parent_id = p_v.get_id();
	}
	else if (current_class != v[i+1].second) {
	  // go down
	  std::cout << "down " << i << "\n";
	  current_parent_class = current_class;
	  current_parent_id = i;	  
	}

	//prev_class = current_class;
	++counter;
	
	++i;
  }
  
  
  return t;
}  
} // namespace pvst
