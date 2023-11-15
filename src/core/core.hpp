#ifndef CORE_HPP
#define CORE_HPP

#include "./constants.hpp"
#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <map>
#include <vector>
#include <iostream>

namespace core {

/*
 * black edge is default
 * TODO: pick a better default
 * gray edge is a bi-edge
 */
enum color { gray, black };

typedef std::pair<std::size_t, std::size_t> size_t_pair;


enum class task_t {
  call, // call variants
  unset // unset
};

enum class input_format_t {
	file_path,
	params,  // CLI params
	unset // unset
};

// TODO: move to source file
// implement << operator for tasks
std::ostream& operator<<(std::ostream& os, const task_t& t);

/*
 * app config
 */
struct config {
  task_t task;
  std::string input_gfa;

  // general
  unsigned char v; // verbosity

  // references
  std::string references_txt; // the path to the file containing the reference paths
  input_format_t ref_input_format;
  std::vector<std::string> reference_paths;


  // -------------
  // Contructor(s)
  // -------------

  // constructor(s)
  config() :
	task(task_t::unset),
	v(0),
	ref_input_format(input_format_t::unset),
	reference_paths(std::vector<std::string>{})
	{}


  // ---------
  // getter(s)
  // ---------
  std::string get_input_gfa() { return this->input_gfa; }
  std::vector<std::string> const& get_reference_paths() const { return this->reference_paths; }
  std::vector<std::string>* get_reference_ptr() { return &this->reference_paths; }
  const std::string& get_references_txt() const { return this->references_txt; }
  std::size_t verbosity() const { return this->v; } // can we avoid this being a size_t?


  // ---------
  // setter(s)
  // ---------
  // as it os from the user not the handlegraph stuff
  void set_ref_input_format(input_format_t f) { this->ref_input_format = f; }
  void add_reference_path(std::string s) { this->reference_paths.push_back(s); }
  void set_reference_paths(std::vector<std::string>&& v) { this->reference_paths = std::move(v); }
  void set_reference_txt_path(std::string&& s) { this->references_txt = std::move(s); }
  void set_references_txt(std::string s) { this->references_txt = s; }
  void set_verbosity(unsigned char v) { this->v = v; }
  void set_input_gfa(std::string s) { this->input_gfa = s; }
  void set_task(task_t t) { this->task = t; }

  // --------
  // Other(s)
  // --------
  void dbg_print() {
	std::cerr << "CLI parameters: " << std::endl;
	std::cerr << "\t" << "verbosity: " << this->verbosity() << "\n";
	std::cerr << "\t" << "task: " << this->task << std::endl;
	std::cerr << "\t" << "input gfa: " << this->input_gfa << std::endl;
	if (this->ref_input_format == input_format_t::file_path) {
	  std::cerr << "\t" << "Reference paths file: " << this->references_txt << std::endl;
	}

	std::cerr << "\t" << "Reference paths (" << this->reference_paths.size() << "): ";
	for (auto it = this->reference_paths.begin(); it != this->reference_paths.end(); ++it) {
	  std::cout << *it;
	  if (std::next(it) != this->reference_paths.end()) { std::cout << ", "; }
	}

	std::cout << std::endl;
	}

};

} // namespace core
#endif
