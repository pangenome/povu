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
#include "../core/core.hpp"



namespace cli {

#define FILE_ERROR(name)                                                       \
  {                                                                            \
	std::string e = "Error, Failed to open the file" + name; \
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
  args::PositionalList<std::string> pathsList(parser, "paths", "files to commit");
	
  parser.Parse();

  //std::cout << "input_gfa: " << bool{input_gfa} << ", value: " << args::get(input_gfa) << std::endl;

  // input gfa is already a c_Str
  app_config.set_input_gfa(args::get(input_gfa));
  
  // either ref list or path list
  // if ref list is not set, then path list must be set
  // -------------

  if (ref_list && std::begin(pathsList) != std::end(pathsList)) {
	std::cerr << "[cli::call_handler]Error: cannot set both ref_list and path_list" << std::endl;
	std::exit(1);

    // throw an argument error
	//throw std::invalid_argument( "[cli::call_handler]Error: cannot set both ref_list and path_list" );
  }
  else if (ref_list) {
	std::cout << "ref_list: " << bool{ref_list} << ", value: " << args::get(ref_list) << std::endl;

	//foo (args::get(ref_list), app_config.get_reference_ptr());
	
	
  }
  else {
	// assume pathsList is not empty
	for (auto &&path : pathsList)
	{
	  app_config.add_reference_path( path);
	}
	  
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
