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
  bool call_variants_;
  std::string input_gfa;
  std::vector<std::string> reference_paths;


  // -----------
  // Contructor(s)
  // -------------

  // constructor(s)
  config() :
	call_variants_(false), reference_paths(std::vector<std::string>{}) {}

  config(std::size_t number_of_refs) :
	call_variants_(false), reference_paths(std::vector<std::string>{}) {
	this->reference_paths.reserve(number_of_refs);
  }

  // --------
  // Getters
  // --------
  std::string get_input_gfa() { return this->input_gfa; }

  bool call_variants() { return this->call_variants_; }

  // --------
  // Setters
  // --------
  // as it os from the user not the handlegraph stuff
  void add_reference_path(std::string s) { this->reference_paths.push_back(s); }

  void set_input_gfa(std::string s) { this->input_gfa = s; }

  void set_call_variants(bool v) { this->call_variants_ = v; }

};

} // namespace core
#endif
