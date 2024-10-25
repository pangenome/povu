#include <cstddef>
#include <cstdlib>
// #include <format>
#include <functional>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <string>
#include <unordered_map>

#include <args.hxx>

#include "./cli.hpp"
#include "app.hpp"

namespace cli {
/**
 * @brief Get the size of a file
 * @param fp file path
 * @return size of the file in bytes
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
  std::size_t file_size = get_file_size(fp);

  v->reserve(file_size);
  read_lines_to_vector_str(fp, v);
  v->shrink_to_fit();
}


void call_handler(args::Subparser &parser, core::config& app_config) {
  args::Group arguments("arguments");
  args::ValueFlag<std::string> input_gfa(parser, "gfa", "path to input gfa [required]", {'i', "input-gfa"}, args::Options::Required);
  args::ValueFlag<std::string> forest_dir(parser, "forest_dir", "dir containing flubble forest [default: .]", {'f', "forest-dir"});
  args::ValueFlag<std::string> output_dir(parser, "output_dir", "Output directory [default: .]", {'o', "output-dir"});
  args::ValueFlag<std::string> ref_list(parser, "ref_list", "path to txt file containing reference haplotypes [optional]", {'p', "path-list"});
  args::ValueFlag<std::string> chrom(parser, "chrom", "graph identifier, default is from GFA file. Chrom column in VCF [optional]", {'c', "chrom"});
  args::Flag undefined_vcf(parser, "undefined_vcf", "Generate VCF file for flubbles without a reference path [default: false]", {'u', "undefined"});
  args::PositionalList<std::string> pathsList(parser, "paths", "list of paths to use as reference haplotypes [optional]");

  parser.Parse();

  app_config.set_task(core::task_t::call);

  // input gfa is already a c_str
  app_config.set_input_gfa(args::get(input_gfa));

  if (forest_dir) {
    app_config.set_forest_dir(args::get(forest_dir));
  }

  if (output_dir) {
    app_config.set_output_dir(args::get(output_dir));
  }

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
    std::cerr << "[cli::call_handler] Error: cannot set both ref_list and path_list" << std::endl;
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


void deconstruct_handler(args::Subparser &parser, core::config& app_config) {
  args::Group arguments("arguments");
  args::ValueFlag<std::string> input_gfa(parser, "gfa", "path to input gfa [required]", {'i', "input-gfa"}, args::Options::Required);
  args::ValueFlag<std::string> output_dir(parser, "output_dir", "Output directory [default: .]", {'o', "output-dir"});

  parser.Parse();
  app_config.set_task(core::task_t::deconstruct);
  // input gfa is already a c_str
  app_config.set_input_gfa(args::get(input_gfa));

  if (output_dir) {
    app_config.set_output_dir(args::get(output_dir));
  }
}


void info_handler(args::Subparser &parser, core::config& app_config) {
  args::Group arguments("arguments");
  args::ValueFlag<std::string> input_gfa(parser, "gfa", "path to input gfa [required]", {'i', "input-gfa"}, args::Options::Required);

  parser.Parse();
  app_config.set_task(core::task_t::info);
  // input gfa is already a c_str
  app_config.set_input_gfa(args::get(input_gfa));
}


int cli(int argc, char **argv, core::config& app_config) {

  args::ArgumentParser p("Use cycle equivalence to call variants");
  args::Group commands(p, "commands");
  args::Command call(commands, "call", "call",
                       [&](args::Subparser &parser) { call_handler(parser, app_config); });
  args::Command deconstruct(commands, "deconstruct", "Find flubbles in the variation graph",
                       [&](args::Subparser &parser) { deconstruct_handler(parser, app_config); });
  args::Command info(commands, "info", "Print information about the graph [use 1 thread for meaningful results]",
                       [&](args::Subparser &parser) { info_handler(parser, app_config); });


  args::Group arguments(p, "arguments", args::Group::Validators::DontCare, args::Options::Global);
  args::Flag version(arguments, "version", "The current version of povu", {"version"});
  args::ValueFlag<int> verbosity(arguments, "verbosity", "Level of output [default: 0]", {'v', "verbosity"});
  args::ValueFlag<int> thread_count(arguments, "threads", "Number of threads to use [default: 1]", {'t', "threads"});
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

  //if (pvst_path) {
  //  app_config.set_pvst_path(args::get(pvst_path));
  //}

  if (args::get(verbosity)) {
    app_config.set_verbosity(args::get(verbosity));
  }

  if (thread_count) {
    app_config.set_thread_count(args::get(thread_count));
  }

  //if (no_sort) {
  //  app_config.set_sort(false);
  //}

  return 0;
}


} // namespace cli
