#ifndef UTILS_HPP
#define UTILS_HPP

#include <chrono>
#include <format>
#include <iostream>
#include <string>
#include <vector>
#include <unordered_map>

namespace povu::utils {

std::string reverse_complement(const std::string& sequence);

/**
 * @brief Concatenates a vector of strings with a given character
 */
// template <typename T> std::string concat_with(const T& v, char delim) {
//   if (v.empty()) { return ""; }

//   std::string s {};
//   for (auto x: v) {
//     s = s + x + delim;
//   }
//   return s.substr(0, s.length()-1);
// }

template <typename Container> std::string concat_with(const Container &v, char delim) {
  std::ostringstream oss;
  auto it = v.begin();
  if (it != v.end()) {
    oss << *it;
    ++it;
    for (; it != v.end(); ++it) {
      oss << delim << *it;
    }
  }
  return oss.str();
}

// TODO rename to print_with_delim or print_with
template <typename T> void print_with_comma(std::ostream& os, const T& v, char delim) {
  if (v.empty()) { return; }

  for (auto it {v.begin()}; it != v.end(); ++it) {
    os << *it;
    if (std::next(it) != v.end()) { os << delim << " "; }
  }
}

/**
  * @brief
 */
void report_time(std::ostream& os, std::string fn_name, std::string action, std::chrono::duration<double> period);


/**
 * @brief Returns the current date in the format YYYYMMDD
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

// TODO : move to povu::types
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

  // get all keys
  std::vector<Key> get_keys() const {
    std::vector<Key> keys;
    for (const auto& [key, _]: keyToValueMap) {
      keys.push_back(key);
    }
    return keys;
  }

  // find key
  bool has_key(const Key& key) const {
    return keyToValueMap.find(key) != keyToValueMap.end();
  }
};

template <typename T> void push_front(std::vector<T>& v, const T& elem) {
  v.insert(v.begin(), elem);
}

/**
 * @brief given a bi-Edged index return its index in the bi-Directed graph
 *
 * @param
 * @param
 * @return
 */
std::size_t to_bidirected_idx(std::size_t be_idx, bool has_dummy=true);

/**
 * @brief given a bi-Directed index return the indexes of the left and right vertices in the bi-Edged graph
 *
 * @param
 * @param
 * @return
 */
std::pair<std::size_t, std::size_t> frm_bidirected_idx(std::size_t bd_idx, bool has_dummy=true);

/**
 * @brief split a string into tokens using a delimiter
 */
void split(const std::string &line, char sep, std::vector<std::string>* tokens);

} // namespace povu::utils
#endif
