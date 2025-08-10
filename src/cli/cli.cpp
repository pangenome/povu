#include <cstddef>
#include <cstdlib>
#include <functional>
#include <fstream>
#include <filesystem>
#include <iostream>
#include <iterator>
#include <string>
#include <unordered_map>

#include <args.hxx> // for command line parsing

#include "./cli.hpp"
#include "app.hpp"

namespace cli {

// TODO: rename ref_list to refs

void call_handler(args::Subparser &parser, core::config& app_config) {
  args::Group arguments("arguments");
  args::ValueFlag<std::string> input_gfa(parser, "gfa", "path to input gfa [required]", {'i', "input-gfa"}, args::Options::Required);
  args::ValueFlag<std::string> forest_dir(parser, "forest_dir", "dir containing flubble forest [default: .]", {'f', "forest-dir"});
  args::ValueFlag<std::string> output_dir(parser, "output_dir", "Output directory [default: .]", {'o', "output-dir"});
  args::ValueFlag<std::string> ref_list(parser, "ref_list", "path to txt file containing reference haplotypes [optional]", {'r', "ref-list"});
  args::ValueFlag<std::string> chrom(parser, "chrom", "graph identifier, default is from GFA file. Chrom column in VCF [optional]", {'c', "chrom"});
  args::ValueFlagList<std::string> path_prefixes(parser, "path_prefix", "All paths beginning with NAME used as reference (multiple allowed) [optional]", {'P', "path-prefix"});
  args::PositionalList<std::string> refsList(parser, "refs", "list of refs to use as reference haplotypes [optional]");

  parser.Parse();

  app_config.set_task(core::task_e::call);

  // input gfa is already a c_str
  app_config.set_input_gfa(args::get(input_gfa));

  if (forest_dir) {
    app_config.set_forest_dir(args::get(forest_dir));
  }

  if (output_dir) {
    app_config.set_output_dir(args::get(output_dir));
  }

  // uses the name of the GFA file
  if (chrom) {
    app_config.set_chrom(std::move(args::get(chrom)));
  }
  else {
    std::filesystem::path filePath(app_config.get_input_gfa());
    app_config.set_chrom(filePath.stem().string());
  }



  /* set graph properties */
  app_config.set_inc_vtx_labels(true);
  app_config.set_inc_refs(true);

  /* One of ref_list, path_prefixes, or list of references must be set—never multiple, and never none */
  int ref_options_set = 0;
  if (ref_list) ref_options_set++;
  if (path_prefixes) ref_options_set++;
  if (std::begin(refsList) != std::end(refsList)) ref_options_set++;

  if (ref_options_set == 0) {
    std::cerr << "[cli::call_handler] Error: need one of: ref_list, path_prefix, or positional refs" << std::endl;
    std::exit(1);
  }
  else if (ref_options_set > 1) {
    std::cerr << "[cli::call_handler] Error: cannot set multiple reference options (ref_list, path_prefix, positional refs)" << std::endl;
    std::exit(1);
  }
  else if (ref_list) {
    app_config.set_ref_input_format(core::input_format_e::file_path);
    app_config.set_reference_txt_path(std::move(args::get(ref_list)));
  }
  else if (path_prefixes) {
    app_config.set_ref_input_format(core::input_format_e::params);
    for (auto &&prefix : args::get(path_prefixes)) {
      app_config.add_path_prefix(prefix);
    }
  }
  else { // we already made sure that the list of refs is not empty
    app_config.set_ref_input_format(core::input_format_e::params);
    // TODO: set this in the call subcommand
    for (auto &&path : refsList) {
      app_config.add_reference_path(path);
    }
  }
}


void deconstruct_handler(args::Subparser &parser, core::config& app_config) {
  args::Group arguments("arguments");
  args::Flag hairpins(parser, "hairpins", "Find hairpins in the variation graph", {'h', "hairpins"});
  args::Flag hubbles(parser, "hubbles", "Find hubbles in the variation graph", {'s', "hubbles"});
  args::ValueFlag<std::string> input_gfa(parser, "gfa", "path to input gfa [required]", {'i', "input-gfa"}, args::Options::Required);
  args::ValueFlag<std::string> output_dir(parser, "output_dir", "Output directory [default: .]", {'o', "output-dir"});

  parser.Parse();
  app_config.set_task(core::task_e::deconstruct);

  if (hairpins) {
    app_config.set_hairpins(true);
  }

  if (hubbles) {
    app_config.set_hubbles(true);
  }

  // input gfa is already a c_str
  app_config.set_input_gfa(args::get(input_gfa));

  if (output_dir) {
    app_config.set_output_dir(args::get(output_dir));
  }
}


void gfa2vcf_handler(args::Subparser &parser, core::config& app_config) {
  args::Group arguments("arguments");
  args::ValueFlag<std::string> input_gfa(parser, "gfa", "path to input gfa [required]", {'i', "input-gfa"}, args::Options::Required);
  args::ValueFlag<std::string> ref_list(parser, "ref_list", "path to txt file containing reference haplotypes [optional]", {'r', "ref-list"});
  args::ValueFlag<std::string> chrom(parser, "chrom", "graph identifier, default is from GFA file. Chrom column in VCF [optional]", {'c', "chrom"});
  args::Flag hairpins(parser, "hairpins", "Find hairpins in the variation graph", {'h', "hairpins"});
  args::Flag hubbles(parser, "hubbles", "Find hubbles in the variation graph", {'s', "hubbles"});
  args::ValueFlagList<std::string> path_prefixes(parser, "path_prefix", "All paths beginning with NAME used as reference (multiple allowed) [optional]", {'P', "path-prefix"});
  args::PositionalList<std::string> refsList(parser, "refs", "list of refs to use as reference haplotypes [optional]");

  parser.Parse();

  app_config.set_task(core::task_e::gfa2vcf);
  app_config.set_input_gfa(args::get(input_gfa));

  if (hairpins) {
    app_config.set_hairpins(true);
  }

  if (hubbles) {
    app_config.set_hubbles(true);
  }

  if (chrom) {
    app_config.set_chrom(std::move(args::get(chrom)));
  }
  else {
    std::filesystem::path filePath(app_config.get_input_gfa());
    app_config.set_chrom(filePath.stem().string());
  }


  app_config.set_inc_vtx_labels(true);
  app_config.set_inc_refs(true);

  int ref_options_set = 0;
  if (ref_list) ref_options_set++;
  if (path_prefixes) ref_options_set++;
  if (std::begin(refsList) != std::end(refsList)) ref_options_set++;

  if (ref_options_set == 0) {
    std::cerr << "[cli::gfa2vcf_handler] Error: need one of: ref_list, path_prefix, or positional refs" << std::endl;
    std::exit(1);
  }
  else if (ref_options_set > 1) {
    std::cerr << "[cli::gfa2vcf_handler] Error: cannot set multiple reference options (ref_list, path_prefix, positional refs)" << std::endl;
    std::exit(1);
  }
  else if (ref_list) {
    app_config.set_ref_input_format(core::input_format_e::file_path);
    app_config.set_reference_txt_path(std::move(args::get(ref_list)));
  }
  else if (path_prefixes) {
    app_config.set_ref_input_format(core::input_format_e::params);
    for (auto &&prefix : args::get(path_prefixes)) {
      app_config.add_path_prefix(prefix);
    }
  }
  else {
    app_config.set_ref_input_format(core::input_format_e::params);
    for (auto &&path : refsList) {
      app_config.add_reference_path(path);
    }
  }
}


void info_handler(args::Subparser &parser, core::config& app_config) {
  args::Group arguments("arguments");
  args::ValueFlag<std::string> input_gfa(parser, "gfa", "path to input gfa [required]", {'i', "input-gfa"}, args::Options::Required);
  args::Flag tips(parser, "tips", "print the tips", {'t', "print_tips"});


  parser.Parse();
  app_config.set_task(core::task_e::info);

  if (tips) {
    app_config.set_print_tips(true);
  }

  // input gfa is already a c_str
  app_config.set_input_gfa(args::get(input_gfa));
}


int cli(int argc, char **argv, core::config& app_config) {

  args::ArgumentParser p("Explore genomic variation in a variation graph");
  args::Group commands(p, "commands");

  args::Command gfa2vcf(commands, "gfa2vcf", "Convert GFA to VCF (decompose + call)",
                       [&](args::Subparser &parser) { gfa2vcf_handler(parser, app_config); });
  args::Command decompose(commands, "decompose", "Find regions of variation",
                       [&](args::Subparser &parser) { deconstruct_handler(parser, app_config); });
  args::Command call(commands, "call", "Generate a VCF from regions of variation",
                       [&](args::Subparser &parser) { call_handler(parser, app_config); });
  args::Command info(commands, "info", "Print graph information [use 1 thread for meaningful results]",
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
