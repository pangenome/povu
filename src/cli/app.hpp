#ifndef CORE_HPP
#define CORE_HPP

#include <cstdint>
#include <string>
#include <unistd.h>
#include <vector>
#include <utility>
#include <iostream>
#include <filesystem>

#include "../../include/common/utils.hpp"

namespace core {
namespace pu = povu::utils;


/*
 * =========
 * App types
 * =========
 */

// rename to task_e
enum class task_e {
  call,        // call variants
  deconstruct, // deconstruct a graph
  info,        // print graph information
  unset        // unset
};

// rename to input_format_e
enum class input_format_e {
  file_path,
  params,  // CLI params
  unset // unset
};
std::ostream& operator<<(std::ostream& os, const task_e& t);

enum class subgraph_category {
  bubble,
  component,
  unset
};

/**
 * @brief app config
 */
struct config {
  task_e task;

  bool inc_hairpins_;
  bool find_hubbles_;
  std::string input_gfa;
  std::filesystem::path forest_dir; //directory containing the flb files
  std::string chrom; // TODO: remove or use
  //std::optional<std::filesystem::path> pvst_path;
  std::filesystem::path output_dir; // output directory for task and deconstruct

  // info
  bool print_tips_; // whether to print tips

  // graph
  bool inc_vtx_labels_; // whether to include vertex labels
  bool inc_refs_; // whether to include references/paths

  // general
  unsigned char v; // verbosity
  bool print_dot_ { true }; // generate dot format graphs

  unsigned int thread_count_ {1}; // number of threads to use

  // references
  std::string references_txt; // the path to the file containing the reference paths
  input_format_e ref_input_format;
  std::vector<std::string> reference_paths; // or just references
  std::vector<std::string> path_prefixes; // path prefixes for reference selection
  bool undefined_vcf; // TODO: remove or use

  // -------------
  // Contructor(s)
  // -------------

  config()
      : task(task_e::unset),
        inc_hairpins_(false),
        find_hubbles_(false),
        forest_dir("."),
        chrom(""), // default is empty string
        output_dir("."), // default is current directory
        print_tips_(false),
        inc_vtx_labels_(false),
        inc_refs_(false),
        v(0),
        thread_count_(1),
        references_txt(""),
        ref_input_format(input_format_e::unset),
        reference_paths(std::vector<std::string>{}),
        path_prefixes(std::vector<std::string>{}),
        undefined_vcf(false)
    {}

  // ---------
  // getter(s)
  // ---------
  bool inc_hairpins() const { return this->inc_hairpins_; }
  bool find_hubbles() const { return this->find_hubbles_; }
  std::string get_input_gfa() const { return this->input_gfa; }
  std::filesystem::path get_forest_dir() const { return this->forest_dir; }
  std::filesystem::path get_output_dir() const { return this->output_dir; }
  bool print_tips() const { return this->print_tips_; }
  bool inc_vtx_labels() const { return this->inc_vtx_labels_; }
  bool inc_refs() const { return this->inc_refs_; }
  const std::string& get_chrom() const { return this->chrom; }
  std::vector<std::string> const& get_reference_paths() const { return this->reference_paths; }
  std::vector<std::string> const& get_path_prefixes() const { return this->path_prefixes; }
  input_format_e get_refs_input_fmt() const { return this->ref_input_format; }
  std::vector<std::string>* get_reference_ptr() { return &this->reference_paths; }
  const std::string& get_references_txt() const { return this->references_txt; }
  std::size_t verbosity() const { return this->v; } // can we avoid this being a size_t?
  unsigned int thread_count() const { return this->thread_count_; }
  bool print_dot() const { return this->print_dot_; }
  bool gen_undefined_vcf() const { return this->undefined_vcf; }
  task_e get_task() const { return this->task; }

  // ---------
  // setter(s)
  // ---------
  // as it os from the user not the handlegraph stuff
  void set_hairpins(bool b) { this->inc_hairpins_ = b; }
  void set_hubbles(bool b) { this->find_hubbles_ = b; }
  void set_chrom(std::string&& s) { this->chrom = s; }
  void set_print_tips(bool b) { this->print_tips_ = b; }
  void set_inc_vtx_labels(bool b) { this->inc_vtx_labels_ = b; }
  void set_inc_refs(bool b) { this->inc_refs_ = b; }
  void set_ref_input_format(input_format_e f) { this->ref_input_format = f; }
  void add_reference_path(std::string s) { this->reference_paths.push_back(s); }
  void add_path_prefix(std::string s) { this->path_prefixes.push_back(s); }
  void set_reference_paths(std::vector<std::string>&& v) { this->reference_paths = std::move(v); }
  void set_path_prefixes(std::vector<std::string>&& v) { this->path_prefixes = std::move(v); }
  void set_reference_txt_path(std::string&& s) { this->references_txt = std::move(s); }
  void set_references_txt(std::string s) { this->references_txt = s; }
  void set_verbosity(unsigned char v) { this->v = v; }
  void set_thread_count(uint8_t t) { this->thread_count_ = t; }
  void set_print_dot(bool b) { this->print_dot_ = b; }
  void set_input_gfa(std::string s) { this->input_gfa = s; }
  void set_forest_dir(std::string s) { this->forest_dir = s; }
  void set_output_dir(std::string s) { this->output_dir = s; }
  void set_task(task_e t) { this->task = t; }
  void set_undefined_vcf(bool b) { this->undefined_vcf = b; }

  // --------
  // other(s)
  // --------
  void dbg_print() {

    const std::string spc = "  "; // could use \t

#ifdef DEBUG
    std::cerr << "povu is in debug mode" << std::endl;
#endif
    std::cerr << "CLI parameters: " << std::endl;

    /* common for all tasks */
    std::cerr << spc << "task: " << this->task << std::endl;
    std::cerr << spc << "verbosity: " << this->verbosity() << "\n";
    std::cerr << spc << "thread count: " << this->thread_count() << "\n";
    std::cerr << spc << "input gfa: " << this->input_gfa << std::endl;
    std::cerr << spc << "output dir: " << this->output_dir << std::endl;

    //std::cerr << spc << "chrom: " << this->chrom << std::endl;
    //std::cerr << spc << "Generate undefined vcf: " << std::boolalpha << this->undefined_vcf << std::endl;

    if (this->get_task() == task_e::call) {
      std::cerr << spc << "forest dir: " << this->forest_dir << std::endl;
      if (this->ref_input_format == input_format_e::file_path) {
        std::cerr << spc << "Reference paths file: " << this->references_txt << std::endl;
      }

      if (!this->path_prefixes.empty()) {
        std::cerr << spc << "Path prefixes (" << this->path_prefixes.size() << "): ";
        pu::print_with_comma(std::cerr, this->path_prefixes, ',');
        std::cerr << std::endl;
      }
      
      std::cerr << spc << "Reference paths (" << this->reference_paths.size() << "): ";
      pu::print_with_comma(std::cerr, this->reference_paths, ',');
      std::cerr << std::endl;

    }
    else if (this->get_task() == task_e::deconstruct) {
#ifdef DEBUG
      std::cerr << spc << "print dot: " << (this->print_dot() ? "yes" : "no") << "\n";
#endif
      std::cerr << spc << "print hairpins: " << (this->inc_hairpins_ ? "yes" : "no") << "\n";
      std::cerr << spc << "find hubbles: " << (this->find_hubbles_ ? "yes" : "no") << "\n";
    }
    else if (this->get_task() == task_e::info) {
      //
    }

    std::cerr << std::endl;
    }
};
} // namespace core
#endif
