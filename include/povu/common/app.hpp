#ifndef CORE_HPP
#define CORE_HPP

#include <cstddef>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <optional>
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

enum class task_e : uint8_t {
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

enum class input_format_e : uint8_t {
	file_path,
	params, // CLI params
	unset	// unset
};
std::ostream &operator<<(std::ostream &os, const task_e &t);

/**
 * @brief app config
 */
struct config {
	task_e task{task_e::unset};

	/* graph decomposition */
	bool inc_hairpins_{false};     // whether to include hairpins
	bool find_subflubbles_{false}; // whether to find subflubbles
	// directory containing the flb files (the forest)
	std::filesystem::path forest_dir{"."};

	/* variation graph config */
	std::string input_gfa{};     // path to input gfa
	bool inc_vtx_labels_{false}; // whether to include vertex labels
	bool inc_refs_{false};	     // whether to include references/paths

	std::size_t chunk_size_{100};
	std::size_t queue_len_{4};

	// general
	unsigned char verbosity_{0}; // verbosity

	// generate dot format graphs
#ifdef DEBUG
	bool print_dot_{true};
	bool print_tips_{false}; // whether to print tips
#else
	bool print_dot_{false};
	bool print_tips_{false}; // whether to print tips
#endif

	unsigned int thread_count_{1}; // number of threads to use

	/* ref handling */

	// the path to the file containing the reference paths
	std::string references_txt{""};
	input_format_e ref_input_format{input_format_e::unset};
	// path prefixes for reference selection
	std::vector<std::string> ref_name_prefixes_{};

	/* VCF file generation */

	// when true, generate a single VCF & write to stdout
	bool stdout_vcf{false};
	// output directory for VCF separate files per ref chosen
	std::filesystem::path output_dir{"."};

	/* genomic region filtering */
	std::optional<std::string> genomic_region_str_{std::nullopt};

	// -------------
	// Contructor(s)
	// -------------

	config() = default;

	// ---------
	// getter(s)
	// ---------
	[[nodiscard]]
	bool inc_hairpins() const
	{
		return this->inc_hairpins_;
	}

	[[nodiscard]]
	bool find_subflubbles() const
	{
		return this->find_subflubbles_;
	}

	[[nodiscard]]
	std::string get_input_gfa() const
	{
		return this->input_gfa;
	}

	[[nodiscard]]
	std::filesystem::path get_forest_dir() const
	{
		return this->forest_dir;
	}

	[[nodiscard]]
	std::filesystem::path get_output_dir() const
	{
		return this->output_dir;
	}

	[[nodiscard]]
	bool print_tips() const
	{
		return this->print_tips_;
	}

	[[nodiscard]]
	bool inc_vtx_labels() const
	{
		return this->inc_vtx_labels_;
	}

	[[nodiscard]]
	bool inc_refs() const
	{
		return this->inc_refs_;
	}

	[[nodiscard]]
	std::size_t get_chunk_size() const
	{
		return this->chunk_size_;
	}

	[[nodiscard]]
	std::size_t get_queue_len() const
	{
		return this->queue_len_;
	}

	[[nodiscard]]
	std::vector<std::string> const &get_ref_name_prefixes() const
	{
		return this->ref_name_prefixes_;
	}

	[[nodiscard]]
	input_format_e get_refs_input_fmt() const
	{
		return this->ref_input_format;
	}

	[[nodiscard]]
	const std::string &get_references_txt() const
	{
		return this->references_txt;
	}

	[[nodiscard]]
	std::size_t verbosity() const
	{
		return this->verbosity_;
	} // can we avoid this being a size_t?

	[[nodiscard]]
	unsigned int thread_count() const
	{
		return this->thread_count_;
	}

	[[nodiscard]]
	bool print_dot() const
	{
		return this->print_dot_;
	}

	[[nodiscard]]
	bool get_stdout_vcf() const
	{
		return this->stdout_vcf;
	}

	[[nodiscard]]
	task_e get_task() const
	{
		return this->task;
	}

	[[nodiscard]]
	const std::optional<std::string> &get_genomic_region() const
	{
		return this->genomic_region_str_;
	}

	[[nodiscard]]
	bool has_genomic_region() const
	{
		return this->genomic_region_str_.has_value();
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

	void add_ref_name_prefix(const std::string &s)
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

	void set_references_txt(const std::string &s)
	{
		this->references_txt = s;
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

	void set_input_gfa(const std::string &s)
	{
		this->input_gfa = s;
	}

	void set_forest_dir(const std::string &s)
	{
		this->forest_dir = s;
	}

	void set_task(task_e t)
	{
		this->task = t;
	}

	void set_output_dir(const std::string &s)
	{
		this->output_dir = s;
	}

	void set_stdout_vcf(bool b)
	{
		this->stdout_vcf = b;
	}

	void set_genomic_region(const std::string &s)
	{
		this->genomic_region_str_ = s;
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
		std::cerr << "CLI parameters:\n";

		/* common for all tasks */
		std::cerr << spc << "task: " << this->task << "\n";
		std::cerr << spc << "verbosity: " << this->verbosity() << "\n";
		std::cerr << spc << "thread count: " << this->thread_count()
			  << "\n";
		std::cerr << spc << "input gfa: " << this->input_gfa << "\n";
		std::cerr << spc << "output dir: " << this->output_dir << "\n";

		if (this->get_task() == task_e::call) {
			std::cerr << spc << "forest dir: " << this->forest_dir
				  << "\n";

			if (this->ref_input_format ==
			    input_format_e::file_path) {
				std::cerr << spc << "Reference paths file: "
					  << this->references_txt << "\n";
			}

			if (!this->ref_name_prefixes_.empty()) {
				std::cerr << spc << "Ref Path prefixes ("
					  << this->ref_name_prefixes_.size()
					  << "): ";
				pu::print_with_comma(std::cerr,
						     this->ref_name_prefixes_,
						     ',');
				std::cerr << "\n";
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

		std::cerr << "\n" << std::flush;
	}
};
} // namespace core
#endif
