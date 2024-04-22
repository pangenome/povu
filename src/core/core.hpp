#ifndef CORE_HPP
#define CORE_HPP


#include <string>
#include <unistd.h>
#include <vector>
#include <utility>
#include <iostream>
#include <optional>
#include <filesystem>


#include "utils.hpp"

namespace core {
/*
 * =============
 * Utility types
 * =============
 */
typedef std::pair<std::size_t, std::size_t> size_t_pair;


/*
 * =========
 * App types
 * =========
 */

// rename to task_e
enum class task_t {
  call, // call variants
  unset // unset
};

// rename to input_format_e
enum class input_format_t {
    file_path,
    params,  // CLI params
    unset // unset
};
std::ostream& operator<<(std::ostream& os, const task_t& t);

/**
 * @brief app config
 */
struct config {
  task_t task;
  std::string input_gfa;
  std::string chrom;
  std::optional<std::filesystem::path> pvst_path;

  // general
  unsigned char v; // verbosity
  bool print_dot_ {false}; // generate dot format graphs

  //
  bool sort_; // sort the graph (default: true)

  // references
  std::string references_txt; // the path to the file containing the reference paths
  input_format_t ref_input_format;
  std::vector<std::string> reference_paths;
  bool undefined_vcf;

  // -------------
  // Contructor(s)
  // -------------

  // constructor(s)
  config() :
    task(task_t::unset),
    chrom(""), // default is empty string
    pvst_path(std::nullopt),
    v(0),
    sort_(true),
    ref_input_format(input_format_t::unset),
    reference_paths(std::vector<std::string>{}),
    undefined_vcf(false)
    {}

  // ---------
  // getter(s)
  // ---------
  std::string get_input_gfa() const { return this->input_gfa; }
  std::optional<std::filesystem::path> get_pvst_path() { return this->pvst_path; }
  const std::string& get_chrom() const { return this->chrom; }
  std::vector<std::string> const& get_reference_paths() const { return this->reference_paths; }
  std::vector<std::string>* get_reference_ptr() { return &this->reference_paths; }
  const std::string& get_references_txt() const { return this->references_txt; }
  std::size_t verbosity() const { return this->v; } // can we avoid this being a size_t?
  bool print_dot() const { return this->print_dot_; }
  bool sort() const { return this->sort_; }
  bool gen_undefined_vcf() const { return this->undefined_vcf; }

  // ---------
  // setter(s)
  // ---------
  // as it os from the user not the handlegraph stuff
  void set_chrom(std::string&& s) { this->chrom = s; }
  void set_ref_input_format(input_format_t f) { this->ref_input_format = f; }
  void add_reference_path(std::string s) { this->reference_paths.push_back(s); }
  void set_reference_paths(std::vector<std::string>&& v) { this->reference_paths = std::move(v); }
  void set_reference_txt_path(std::string&& s) { this->references_txt = std::move(s); }
  void set_references_txt(std::string s) { this->references_txt = s; }
  void set_verbosity(unsigned char v) { this->v = v; }
  void set_sort(bool b) { this->sort_ = b; }
  void set_print_dot(bool b) { this->print_dot_ = b; }
  void set_input_gfa(std::string s) { this->input_gfa = s; }
  void set_pvst_path(std::string s) { this->pvst_path = s; }
  void set_task(task_t t) { this->task = t; }
  void set_undefined_vcf(bool b) { this->undefined_vcf = b; }

  // --------
  // Other(s)
  // --------
  void dbg_print() {
    std::cerr << "CLI parameters: " << std::endl;
    std::cerr << "\t" << "verbosity: " << this->verbosity() << "\n";
    std::cerr << "\t" << "print dot: " << (this->print_dot() ? "yes" : "no") << "\n";
    std::cerr << "\t" << "sort: " << (this->sort() ? "yes" : "no") << "\n";
    std::cerr << "\t" << "task: " << this->task << std::endl;
    std::cerr << "\t" << "input gfa: " << this->input_gfa << std::endl;
    std::cerr << "\t" << "chrom: " << this->chrom << std::endl;
    std::cerr << "\t" << "Generate undefined vcf: " << std::boolalpha << this->undefined_vcf << std::endl;
    if (this->ref_input_format == input_format_t::file_path) {
      std::cerr << "\t" << "Reference paths file: " << this->references_txt << std::endl;
    }

    std::cerr << "\t" << "Reference paths (" << this->reference_paths.size() << "): ";
    utils::print_with_comma(std::cerr, this->reference_paths, ',');


    std::cerr << std::endl;
    }
};


} // namespace core
#endif
