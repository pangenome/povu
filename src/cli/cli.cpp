#include <chrono>
#include <iostream>
#include <iterator>
#include <string>
#include <unordered_map>
#include <functional>
#include <args.hxx>

#include "../core/core.hpp"

namespace cli {

// version string constant
const std::string VERSION = "0.0.0-alpha";
  

  void call_handler(args::Subparser &parser, core::config& app_config) {
  args::Group arguments("arguments");
  args::ValueFlag<std::string> input_gfa(parser, "gfa", "path to input gfa [required]", {'i', "input-gfa"}, args::Options::Required);
  args::ValueFlag<std::string> ref_list(parser, "ref_list", "path to txt file containing reference haplotypes [optional]", {'p', "path-list"});
  args::PositionalList<std::string> pathsList(parser, "paths", "files to commit");
	
  parser.Parse();

  std::cout << "input_gfa: " << bool{input_gfa} << ", value: " << args::get(input_gfa) << std::endl;

  // either ref list or path list
  // if ref list is not set, then path list must be set
  // -------------
		
  std::cout << "ref_list: " << bool{ref_list} << ", value: " << args::get(ref_list) << std::endl;
	
  for (auto &&path : pathsList)
  {
	std::cout << ' ' << path;
  }

  std::cout << std::endl;

  //if (message)
  //{
  //    std::cout << "message: " << args::get(message) << std::endl;
  //}
}


  
int cli(int argc, char **argv, core::config& app_config) {

  args::ArgumentParser p("Use cycle equivalence call variants");
  args::Group commands(p, "commands");
  args::Command commit(commands, "call", "call",
					   [&](args::Subparser &parser) { call_handler(parser, app_config); });
  //p.Prog(argv[0]);
  args::Group arguments(p, "arguments", args::Group::Validators::DontCare, args::Options::Global);
  args::Flag version(arguments, "version", "show the ...", {"version"});
  args::ValueFlag<int> verbosity(arguments, "verbosity", "Level out output", {'v', "verbosity"});
  args::HelpFlag h(arguments, "help", "help", {'h', "help"});


  try
  {
	p.ParseCLI(argc, argv);
  }
  catch (args::Help)
  {
	std::cout << p;
  }
  catch (args::Error& e)
  {
	std::cerr << e.what() << std::endl << p;
	return 1;
  }

  if (verbosity) { std::cout << "i: " << args::get(verbosity) << std::endl; }
  if (version) { std::cout << VERSION << std::endl; }


  return 0;
}


} // namespace cli
