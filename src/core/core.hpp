#ifndef CORE_HPP
#define CORE_HPP

#include "./constants.hpp"
#include <string>
#include <utility>
#include <map>
#include <vector>


namespace core {
  
/*
 * black edge is default
 * gray edge is a bi-edge
 */
enum color { gray, black };

typedef std::pair<std::size_t, std::size_t> size_t_pair;

// app config

struct config {
  std::vector<std::string> reference_paths;

  // contructor
  config():  reference_paths(std::vector<std::string>{}) {
  }
  
  config(std::size_t number_of_refs) :  reference_paths(std::vector<std::string>{}) {
	this->reference_paths.reserve(number_of_refs);
  }

  // as it os from the user not the handlegraph stuff
  void add_reference_path(std::string s) { this->reference_paths.push_back(s); }
  
};
  
} // namespace core
#endif
