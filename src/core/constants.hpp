#ifndef CONSTANTS_HPP
#define CONSTANTS_HPP

#include <string>
#include <limits>



namespace core::constants {

// colors
const std::string gray{"gray"};
const std::string black{"black"};
const std::string red{"red"};

// numeric
const std::size_t SIZE_T_MIN = std::numeric_limits<size_t>::min();
const std::size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();
const std::size_t UNDEFINED_SIZE_T = std::numeric_limits<size_t>::max();
const int UNDEFINED_INT = std::numeric_limits<int>::min();

// strings
const std::string EMPTY_SET = "\u2205";
const std::string UNDEFINED_VALUE = "\u2205";
}
#endif
