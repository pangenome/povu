#include <iostream>

#include "./utils.hpp"

namespace utils {

void print_with_comma(std::unordered_set<std::size_t>& iterable) {
  for (auto it = iterable.begin(); it != iterable.end(); ++it) {
	std::cerr << *it;

	if (std::next(it) != iterable.end()){ std::cerr << ", "; }
  }
}

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
} // namespace utils
