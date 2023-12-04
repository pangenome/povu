#include <iostream>
#include <string>
#include <numeric>
#include <ctime>
#include <iomanip>

#include "./utils.hpp"

namespace utils {

void print_with_comma(std::unordered_set<std::size_t>& iterable) {
  for (auto it = iterable.begin(); it != iterable.end(); ++it) {
	std::cerr << *it;

	if (std::next(it) != iterable.end()){ std::cerr << ", "; }
  }
}

std::string concat_with (const std::vector<std::string>& v, char c) {
  return std::accumulate(v.begin(), v.end(), std::string(),
						 [&c](const std::string& a, const std::string& b) -> std::string {
						   return a + c + b;
						 }
	);
};

char complement(char nucleotide) {
  switch (nucleotide) {
  case 'A': return 'T';
  case 'T': return 'A';
  case 'C': return 'G';
  case 'G': return 'C';
	// You might want to handle cases like 'N' for unknown nucleotides
  default: return nucleotide;
  }
}

std::string reverse_complement(const std::string& sequence) {
  std::string rc_sequence;

  for (auto it = sequence.rbegin(); it != sequence.rend(); ++it) {
	rc_sequence += complement(*it);
  }

  return rc_sequence;
}

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

std::vector<std::string> immutable_erase(std::vector<std::string>& v, std::size_t idx) {
  std::vector<std::string> v_ = v;
  v_.erase(v_.begin()+idx, v_.begin()+idx);
  return v_;
}
  
} // namespace utils
