#include <iostream>

#include "./utils.hpp"

namespace utils {
  void print_with_comma(std::unordered_set<std::size_t>& iterable) {
  for (auto it = iterable.begin(); it != iterable.end(); ++it) {
	std::cerr << *it;

	if (std::next(it) != iterable.end()){ std::cerr << ", "; }
  }
}
} // namespace utils
