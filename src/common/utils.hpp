#ifndef UTILS_HPP
#define UTILS_HPP

#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>
#include <unordered_set>

namespace utils {

// TODO: - generalize for other iterators
//       - pass os stream to print to
void print_with_comma(std::unordered_set<std::size_t>& iterable);
void print_with_comma(std::unordered_set<id_t>&& iterable);
std::string reverse_complement(const std::string& sequence);

/**
 * Concatenates a vector of strings with a given character
 */
//std::string concat_with (const std::vector<std::string>& v, char c);
//template <typename T> std::string concat_with(const T& v, char delim);

template <typename T> std::string concat_with(const T& v, char delim) {
  if (v.empty()) { return ""; }

  std::string s {};
  for (auto x: v) { s = s + x + delim; }
  return s.substr(0, s.length()-1);
}

template <typename T> void print_with_comma(std::ostream& os, const T& v, char delim) {
  if (v.empty()) { return; }

  for (auto it {v.begin()}; it != v.end(); ++it) {
    os << *it;
    if (std::next(it) != v.end()) { os << delim << " "; }
  }
}

/**
 * Returns the current date in the format YYYYMMDD
 */
std::string today();

/**
 * @brief
 *
 * @param v: the vector whose value is to be erased passed by copy to avoid mutating the original
 * @param idx: the index to be erased
 * @return a vector with the value at the given index erased
 */
std::vector<std::string> immutable_erase(std::vector<std::string> v, std::size_t idx);

template <typename Key, typename Value>
class TwoWayMap {
    std::unordered_map<Key, Value> keyToValueMap;
    std::unordered_map<Value, Key> valueToKeyMap;

public:
    // Insert a key-value pair
    void insert(const Key& key, const Value& value) {
        keyToValueMap[key] = value;
        valueToKeyMap[value] = key;
    }

    // Lookup value by key
    Value get_value(const Key& key) const {
        auto it = keyToValueMap.find(key);
        if (it != keyToValueMap.end()) {
            return it->second;
        }
        // Return a default-constructed Value or throw an exception, depending on your use case
        return Value{};
    }

    // Lookup key by value
    Key get_key(const Value& value) const {
        auto it = valueToKeyMap.find(value);
        if (it != valueToKeyMap.end()) {
            return it->second;
        }
        // Return a default-constructed Key or throw an exception, depending on your use case
        return Key{};
    }
};

struct unordered_pair{
  std::size_t l;
  std::size_t r;

  unordered_pair(std::size_t l,std::size_t r):
    l(std::min(l,r)), r(std::max(l,r)) {}

  // spaceship operator
  friend constexpr auto operator<=>(unordered_pair, unordered_pair) = default;
};



} // namespace utils
#endif
