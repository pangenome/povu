#include <chrono>
#include <cstdlib>
#include <iostream>
#include <iterator>
#include <string>
#include <unordered_map>
#include <functional>
#include <stdexcept>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <args.hxx>
#include <filesystem>



#include "./cli.hpp"

namespace cli {

#define FILE_ERROR(name)                                                       \
  {                                                                            \
    std::string e = "Error, Failed to open the file " + name; \
    throw std::invalid_argument(e);                                              \
  }



// version string constant
const std::string VERSION = "0.0.0-alpha";

/**
 * Get the size of a file
 */
std::size_t get_file_size(const std::string& fp) {
  std::streampos begin, end;
  std::ifstream f (fp);

  if (!f) { FILE_ERROR(fp); }

  begin = f.tellg();
  f.seekg (0, std::ios::end);
  end = f.tellg();
  f.close();

  return end-begin;
}


/**
 * read the entire file into a string
 * should be faster for small inputs to read the entire file into a string and
 * process it at once
 * will perform whitespace normalization/formatting
 */
void read_lines_to_vector_str(const std::string& fp, std::vector<std::string> *v) {
  std::ifstream f{fp};
  std::string temp;

  if (!f) { FILE_ERROR(fp); }

  while (f >> temp) { v->push_back(temp); }
}

void fp_to_vector (const std::string& fp, std::vector<std::string>* v) {
  size_t file_size = get_file_size(fp);

  v->reserve(file_size);
  read_lines_to_vector_str(fp, v);
  v->shrink_to_fit();
}

void call_handler(args::Subparser &parser, core::config& app_config) {
  args::Group arguments("arguments");
  args::ValueFlag<std::string> input_gfa(parser, "gfa", "path to input gfa [required]", {'i', "input-gfa"}, args::Options::Required);
  args::ValueFlag<std::string> ref_list(parser, "ref_list", "path to txt file containing reference haplotypes [optional]", {'p', "path-list"});
  args::ValueFlag<std::string> chrom(parser, "chrom", "graph identifier, default is from GFA file. Chrom column in VCF [optional]", {'c', "chrom"});
  args::Flag undefined_vcf(parser, "undefined_vcf", "Generate VCF file for flubbles without a reference path [default: false]", {'u', "undefined"});
  args::PositionalList<std::string> pathsList(parser, "paths", "list of paths to use as reference haplotypes [optional]");

  parser.Parse();

  app_config.set_task(core::task_t::call);

  // input gfa is already a c_str
  app_config.set_input_gfa(args::get(input_gfa));

  if (chrom) {
    app_config.set_chrom(std::move(args::get(chrom)));
  }
  else {
    std::filesystem::path filePath(app_config.get_input_gfa());
    app_config.set_chrom(filePath.stem().string());
  }

  if (undefined_vcf) {
    app_config.set_undefined_vcf(true);
  }

  // either ref list or path list
  // if ref list is not set, then path list must be set
  // -------------

  if (ref_list && std::begin(pathsList) != std::end(pathsList)) {
    std::cerr << "[cli::call_handler]Error: cannot set both ref_list and path_list" << std::endl;
    std::exit(1);
  }
  else if (ref_list) {
    app_config.set_ref_input_format(core::input_format_t::file_path);
    app_config.set_reference_txt_path(std::move(args::get(ref_list)));
    std::vector<std::string> ref_paths;
    read_lines_to_vector_str(app_config.get_references_txt(), &ref_paths);
    app_config.set_reference_paths(std::move(ref_paths));
  }
  else {
    app_config.set_ref_input_format(core::input_format_t::params);
    // assume pathsList is not empty
    for (auto &&path : pathsList) {
      app_config.add_reference_path(path);
    }
  }
}



int cli(int argc, char **argv, core::config& app_config) {

  args::ArgumentParser p("Use cycle equivalence to call variants");
  args::Group commands(p, "commands");
  args::Command commit(commands, "call", "call",
                       [&](args::Subparser &parser) { call_handler(parser, app_config); });

  args::Group arguments(p, "arguments", args::Group::Validators::DontCare, args::Options::Global);
  args::Flag version(arguments, "version", "The current version of povu", {"version"});
  args::ValueFlag<int> verbosity(arguments, "verbosity", "Level of output", {'v', "verbosity"});
  args::ValueFlag<std::string> pvst_path(arguments, "pvst_file_name", "PVST output file path", {'t', "pvst-path"});
  args::Flag no_sort(arguments, "no_sort", "Disable sorting (not recommended)", {'n', "no-sort"});
  args::HelpFlag h(arguments, "help", "help", {'h', "help"});


  try {
    p.ParseCLI(argc, argv);
  }
  catch (args::Help& _) {
    std::cout << p;
  }
  catch (args::Error& e) {
    // only run this if the user is not requesting to print the version
    if (!version) {
      std::cerr << e.what() << std::endl << p;
      return 1;
    }
  }

  if (version) {
    std::cout << VERSION << std::endl;
    std::exit(0);
  }

  if (pvst_path) {
    app_config.set_pvst_path(args::get(pvst_path));
  }

  if (args::get(verbosity)) {
    app_config.set_verbosity(args::get(verbosity));
  }

  if (no_sort) {
    app_config.set_sort(false);
  }

  return 0;
}


} // namespace cli
