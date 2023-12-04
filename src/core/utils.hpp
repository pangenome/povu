#ifndef UTILS_HPP
#define UTILS_HPP

#include <unordered_set>
#include <string>
#include <vector>

namespace utils {
  
// TODO: - generalize for other iterators
//       - pass os stream to print to
void print_with_comma(std::unordered_set<std::size_t>& iterable);
std::string reverse_complement(const std::string& sequence);
/**
 * Concatenates a vector of strings with a given character
 */ 
std::string concat_with (const std::vector<std::string>& v, char c);

/**
 * Returns the current date in the format YYYYMMDD
 */
std::string today();

/**
 * @brief
 *
 * @param v the vector whose value is to be erased
 * @param idx the index to be erased
 * @return a vector with the value at the given index erased
 */
std::vector<std::string> immutable_erase(std::vector<std::string>& v, std::size_t idx);
} // namespace utils
#endif
