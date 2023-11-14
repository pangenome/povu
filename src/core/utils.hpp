#ifndef UTILS_HPP
#define UTILS_HPP

#include <unordered_set>

namespace utils {
  
// TODO: - generalize for other iterators
//       - pass os stream to print to
void print_with_comma(std::unordered_set<std::size_t>& iterable);

} // namespace utils
#endif
