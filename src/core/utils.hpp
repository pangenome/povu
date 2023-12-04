#ifndef UTILS_HPP
#define UTILS_HPP

#include <unordered_set>
#include <string>

namespace utils {
  
// TODO: - generalize for other iterators
//       - pass os stream to print to
void print_with_comma(std::unordered_set<std::size_t>& iterable);
std::string reverse_complement(const std::string& sequence);
} // namespace utils
#endif
