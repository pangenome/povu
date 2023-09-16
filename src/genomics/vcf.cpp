
#include <cstddef>
#include <iostream>
#include <stack>
#include <queue>
#include <string>
#include <map>
#include <set>
#include <unordered_set>
#include <vector>
#include <utility>
#include <ctime>
#include <iomanip>


#include "../pvst/pvst.hpp"
#include "../graph/tree.hpp"
#include "../graph/digraph.hpp"
#include "../core/constants.hpp"
#include "../core/core.hpp"
#include "./genomics.hpp"

namespace vcf {
  
std::string today() {
    // Get the current time
    std::time_t currentTime = std::time(nullptr);

    // Convert the current time to a local time structure
    std::tm* localTime = std::localtime(&currentTime);

    // Extract date and time components
    int year = localTime->tm_year + 1900; // Years since 1900
    int month = localTime->tm_mon + 1;    // Months are 0-based
    int day = localTime->tm_mday;         // Day of the month

    // Create a stringstream to store the formatted date and time
    std::stringstream dateTimeStream;
    dateTimeStream << year
                   << std::setw(2) << std::setfill('0') << month 
                   << std::setw(2) << std::setfill('0') << day;

    return dateTimeStream.str();
}

void print_vcf_header(std::string ref_name) {
std::string blank_str = "xxxx";
  std::string sample_name = blank_str;
  
  std::cout << "##fileformat=VCFv4.2" << std::endl;
  std::cout << "##fileDate=" << today() << std::endl;
  std::cout << "##source=flubble" << std::endl;
  std::cout << "##reference=" << ref_name << std::endl;
  std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + blank_str << std::endl;
}

// path idx, path pos 
std::pair<std::size_t, std::size_t>
get_ref_path_idx(std::vector<std::vector<std::size_t>> paths,
				 digraph::DiGraph dg,
				 std::size_t ref_path_id) {

  if (false) {
	std::cout << "[vcf::get_ref_path_idx]" << std::endl;
  }
  
  // find the reference path
  for (std::size_t i{}; i < paths.size(); ++i) {
	std::vector<std::size_t> path = paths[i];

	// only check the front of the path
	// check the next after the front of the path which is the flubble start
	std::size_t n_id = path.front();

	//std::cout << "n id "<< n_id << std::endl;
	
	for (auto p : dg.get_vertex(n_id).get_paths()) {
	  // std::cout << "path id: " << p.first << " path start "<< p.second<< std::endl;
	  if (p.first == ref_path_id) {

		return std::make_pair(i, p.second);
	  }
	}
  }
  
  // TODO: throw an exception? or default to zero?
  return std::make_pair(
	core::constants::UNDEFINED_SIZE_T, core::constants::UNDEFINED_SIZE_T);
}

/**
 * return code:
 * - < 0 means ...
 * - 0 means success
 * - 1 means that the reference path is not in the flubble
 */
int extract_vcf_variants(std::vector<std::vector<std::size_t>> paths,
						 digraph::DiGraph dg,
						 std::size_t ref_path_id,
						 std::string& ref_path_str,
						 std::vector<std::string>& alt_path_strs,
						 std::size_t& ref_path_pos) {
  
  std::pair<std::size_t, std::size_t> idx_n_pos =
    get_ref_path_idx(paths, dg, ref_path_id);

  if (ref_path_id != core::constants::UNDEFINED_SIZE_T
	  && idx_n_pos.first == core::constants::UNDEFINED_SIZE_T
	  && idx_n_pos.second == core::constants::UNDEFINED_SIZE_T) {	
	return -1;  
  }
  
  ref_path_pos = idx_n_pos.second;


  for (std::size_t i{}; i < paths.size(); ++i) {
	std::vector<std::size_t> path = paths[i];

	// the first and last elements are flubble boundaries
	for (std::size_t j{1}; j < path.size()-1; ++j) {
	  std::size_t n = path[j];

	  if (ref_path_id != core::constants::UNDEFINED_SIZE_T && i == idx_n_pos.first) {
		ref_path_str += dg.get_vertex(n).get_seq();
	  }
	  else {
		alt_path_strs[i] += dg.get_vertex(n).get_seq();
	  }
	}
  }
  
  if (ref_path_id == core::constants::UNDEFINED_SIZE_T) {	
	return 1;  
  } else {
	return 0;	  
  }
}

std::string join_strs (std::vector<std::string> strVector) {
  std::string joinedString;
  for (size_t i = 0; i < strVector.size(); ++i) {
	joinedString += strVector[i];
	if (i > 0 && i < strVector.size() - 1) { joinedString += ","; }
  }
  return joinedString;
}

// prints a single line of the vcf file
// called for each bubble
// if paths in that flubble got printed then return true
bool print_vcf(std::vector<std::vector<std::size_t>> paths, digraph::DiGraph dg,
			   std::string ref_path_name, std::size_t ref_path_id) {

  std::string ref_path_str{};
  std::vector<std::string> alt_path_strs(paths.size(), std::string{});

  std::size_t ref_path_pos{};

  // leave early in case of error
  int res =
	extract_vcf_variants(
	  paths, dg, ref_path_id, ref_path_str, alt_path_strs, ref_path_pos);
	
  if (res < 0) { return false; }

  std::string blank_str = "xxxx";

  std::string chrom = 
	ref_path_name == "undefined" ? "undef" : ref_path_name;
  
  std::string pos =
	ref_path_pos == core::constants::UNDEFINED_SIZE_T ? "-1" : std::to_string( ref_path_pos);
  
  std::string id = blank_str;

  std::string ref = ref_path_str;
  std::string alt = join_strs(alt_path_strs);
  std::string qual = blank_str;
  std::string filter = blank_str;
  std::string info = blank_str;
  std::string format = blank_str;
  
  std::cout << chrom << "\t" << pos << "\t" << id << "\t"
			<< ref << "\t" << alt << "\t"
			<< qual << "\t" << filter << "\t" << info << "\t" << format
			<< "\n";

  return true;
}
} // namespace vcf
