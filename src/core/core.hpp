#ifndef CORE_HPP
#define CORE_HPP

#include "./constants.hpp"
#include <utility>

namespace core {
  
/*
 * black edge is default
 * gray edge is a bi-edge
 */
enum color { gray, black };

typedef std::pair<std::size_t, std::size_t> size_t_pair;
  
} // namespace core
#endif
