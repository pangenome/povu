
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

    // Get the formatted date and time as a string
    //std::string dateTimeString = dateTimeStream.str();

    //std::cout << dateTimeString << "\n";

    return dateTimeStream.str();
}

void print_vcf_header() {
std::string blank_str = "xxxx";
  std::string sample_name = blank_str;
  
  std::cout << "##fileformat=VCFv4.2" << std::endl;
  std::cout << "##fileDate=" << today() << std::endl;
  std::cout << "##source=flubble" << std::endl;
  std::cout << "##reference=" + blank_str << std::endl;
  std::cout << "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + blank_str << std::endl;
  
}
  
void print_vcf_variants(std::vector<std::vector<std::size_t>> paths, digraph::DiGraph dg) {

  std::string blank_str = "xxxx";

  std::string chrom = blank_str;
  std::string pos = blank_str;
  std::string id = blank_str;

  std::string ref = blank_str;
  std::string alt = blank_str;
  std::string qual = blank_str;
  std::string filter = blank_str;
  std::string info = blank_str;
  std::string format = blank_str;
  
  
  
  std::cout << chrom << "\t" << pos << "\t" << id << "\t"; 
  // path 1 is the reference
  for (std::size_t i{}; i < paths.size(); ++i) {
	std::vector<std::size_t> path = paths[i];
	std::string path_str = "";
	for (std::size_t j{1}; j < path.size()-1; ++j) {
	  std::size_t n = path[j];
	  //std::cout << n << "\n";
	  path_str += dg.get_vertex(n).get_seq();
	  //path_str += std::to_string(v);
	}
	
	std::cout << path_str;
	if (i > 0 && i < paths.size()-1) {
	  std::cout << ",";
	}
	else { // i == 0
	  std::cout << "\t";
	}
	
  }

  std::cout << qual << "\t" << filter << "\t" << info << "\t" << format;

  std::cout << "\n";
}
} // namespace vcf
