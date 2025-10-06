#ifndef CORE_HPP
#define CORE_HPP

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <ostream>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

#include "povu/common/log.hpp"
#include "povu/common/utils.hpp"

namespace core
{
constexpr std::string_view MODULE = "povu::app::core";

namespace pu = povu::utils;

/*
 * =========
 * App types
 * =========
 */

// rename to task_e
enum class task_e {
	call,	   // call variants
	decompose, // deconstruct a graph
	gfa2vcf,   // convert GFA directly to VCF
	info,	   // print graph information
	unset	   // unset
};

inline const char *to_str(task_e t)
{
	switch (t) {
	case task_e::call:
		return "call";
	case task_e::decompose:
		return "decompose";
	case task_e::gfa2vcf:
		return "gfa2vcf";
	case task_e::info:
		return "info";
	default:
		return "unset";
	}
}

inline std::ostream &operator<<(std::ostream &os, const task_e &t)
{
	return os << to_str(t);
}

// rename to input_format_e
enum class input_format_e {
	file_path,
	params, // CLI params
	unset	// unset
};
std::ostream &operator<<(std::ostream &os, const task_e &t);

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
	bool find_subflubbles_;
	std::string input_gfa;
	std::filesystem::path forest_dir; // directory containing the flb files
	std::string chrom;		  // TODO: remove or use
	// std::optional<std::filesystem::path> pvst_path;
	std::filesystem::path
		output_dir; // output directory for task and deconstruct

	// info
	bool print_tips_; // whether to print tips

	// graph
	bool inc_vtx_labels_; // whether to include vertex labels
	bool inc_refs_;	      // whether to include references/paths

	std::size_t chunk_size_{100};
	std::size_t queue_len_{4};

	// general
	bool prog_{false};	  // show progress bars
	unsigned char verbosity_; // verbosity
	bool print_dot_{true};	  // generate dot format graphs

	unsigned int thread_count_{1}; // number of threads to use

	// references
	std::string references_txt; // the path to the file containing the
				    // reference paths
	input_format_e ref_input_format;
	std::vector<std::string>
		ref_name_prefixes_; // path prefixes for reference selection
	bool stdout_vcf; // output single VCF to stdout instead of separate
			 // files

	// -------------
	// Contructor(s)
	// -------------

	config()
	    : task(task_e::unset), inc_hairpins_(false),
	      find_subflubbles_(false), forest_dir("."),
	      chrom(""),       // default is empty string
	      output_dir("."), // default is current directory
	      print_tips_(false), inc_vtx_labels_(false), inc_refs_(false),
	      verbosity_(0), thread_count_(1), references_txt(""),
	      ref_input_format(input_format_e::unset),
	      ref_name_prefixes_(std::vector<std::string>{}), stdout_vcf(false)
	{}

	// ---------
	// getter(s)
	// ---------
	bool inc_hairpins() const
	{
		return this->inc_hairpins_;
	}

	bool find_subflubbles() const
	{
		return this->find_subflubbles_;
	}

	std::string get_input_gfa() const
	{
		return this->input_gfa;
	}

	std::filesystem::path get_forest_dir() const
	{
		return this->forest_dir;
	}

	std::filesystem::path get_output_dir() const
	{
		return this->output_dir;
	}

	bool print_tips() const
	{
		return this->print_tips_;
	}

	bool inc_vtx_labels() const
	{
		return this->inc_vtx_labels_;
	}

	bool inc_refs() const
	{
		return this->inc_refs_;
	}

	std::size_t get_chunk_size() const
	{
		return this->chunk_size_;
	}

	std::size_t get_queue_len() const
	{
		return this->queue_len_;
	}

	std::vector<std::string> const &get_ref_name_prefixes() const
	{
		return this->ref_name_prefixes_;
	}

	input_format_e get_refs_input_fmt() const
	{
		return this->ref_input_format;
	}

	const std::string &get_references_txt() const
	{
		return this->references_txt;
	}

	bool show_progress() const
	{
		return this->prog_;
	}

	std::size_t verbosity() const
	{
		return this->verbosity_;
	} // can we avoid this being a size_t?

	unsigned int thread_count() const
	{
		return this->thread_count_;
	}

	bool print_dot() const
	{
		return this->print_dot_;
	}

	bool get_stdout_vcf() const
	{
		return this->stdout_vcf;
	}

	task_e get_task() const
	{
		return this->task;
	}

	// ---------
	// setter(s)
	// ---------
	// as it os from the user not the handlegraph stuff
	void set_hairpins(bool b)
	{
		this->inc_hairpins_ = b;
	}

	void set_subflubbles(bool b)
	{
		this->find_subflubbles_ = b;
	}

	void set_chunk_size(std::size_t s)
	{
		this->chunk_size_ = s;
	}

	void set_queue_len(std::size_t l)
	{
		this->queue_len_ = l;
	}

	void set_print_tips(bool b)
	{
		this->print_tips_ = b;
	}

	void set_inc_vtx_labels(bool b)
	{
		this->inc_vtx_labels_ = b;
	}

	void set_inc_refs(bool b)
	{
		this->inc_refs_ = b;
	}

	void set_ref_input_format(input_format_e f)
	{
		this->ref_input_format = f;
	}

	void add_ref_name_prefix(std::string s)
	{
		this->ref_name_prefixes_.push_back(s);
	}

	void set_ref_name_prefixes(std::vector<std::string> &&v)
	{
		this->ref_name_prefixes_ = std::move(v);
	}

	void set_reference_txt_path(std::string &&s)
	{
		this->references_txt = std::move(s);
	}

	void set_references_txt(std::string s)
	{
		this->references_txt = s;
	}

	void set_progress(bool b)
	{
		this->prog_ = b;
	}

	void set_verbosity(unsigned char v)
	{
		this->verbosity_ = v;
	}

	void set_thread_count(uint8_t t)
	{
		this->thread_count_ = t;
	}

	void set_print_dot(bool b)
	{
		this->print_dot_ = b;
	}

	void set_input_gfa(std::string s)
	{
		this->input_gfa = s;
	}

	void set_forest_dir(std::string s)
	{
		this->forest_dir = s;
	}

	void set_output_dir(std::string s)
	{
		this->output_dir = s;
	}

	void set_task(task_e t)
	{
		this->task = t;
	}

	void set_stdout_vcf(bool b)
	{
		this->stdout_vcf = b;
	}

	// --------
	// other(s)
	// --------
	void dbg_print()
	{

		const std::string spc = "  "; // could use \t

#ifdef DEBUG
		INFO("povu is in debug mode");
#endif
		std::cerr << "CLI parameters: " << std::endl;

		/* common for all tasks */
		std::cerr << spc << "task: " << this->task << std::endl;
		std::cerr << spc << "verbosity: " << this->verbosity() << "\n";
		std::cerr << spc << "thread count: " << this->thread_count()
			  << "\n";
		std::cerr << spc << "input gfa: " << this->input_gfa
			  << std::endl;
		std::cerr << spc << "output dir: " << this->output_dir
			  << std::endl;

		if (this->get_task() == task_e::call) {
			std::cerr << spc << "forest dir: " << this->forest_dir
				  << std::endl;

			if (this->ref_input_format ==
			    input_format_e::file_path) {
				std::cerr << spc << "Reference paths file: "
					  << this->references_txt << std::endl;
			}

			if (!this->ref_name_prefixes_.empty()) {
				std::cerr << spc << "Ref Path prefixes ("
					  << this->ref_name_prefixes_.size()
					  << "): ";
				pu::print_with_comma(std::cerr,
						     this->ref_name_prefixes_,
						     ',');
				std::cerr << std::endl;
			}
		}
		else if (this->get_task() == task_e::decompose) {
#ifdef DEBUG
			std::cerr << spc << "print dot: "
				  << (this->print_dot() ? "yes" : "no") << "\n";
#endif
			std::cerr << spc << "print hairpins: "
				  << (this->inc_hairpins_ ? "yes" : "no")
				  << "\n";
			std::cerr << spc << "find subflubbles: "
				  << (this->find_subflubbles() ? "yes" : "no")
				  << "\n";
		}
		else if (this->get_task() == task_e::info) {
			//
		}

		std::cerr << std::endl;
	}
};
} // namespace core
#endif
