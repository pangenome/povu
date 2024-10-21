#include "./app.hpp"

namespace core {
std::ostream& operator<<(std::ostream& os, const task_t& t)  {
  switch (t) {
    case task_t::call:
      os << "call";
      break;
    case task_t::deconstruct:
      os << "deconstruct";
      break;
    case task_t::info:
      os << "info";
      break;
    default:
      os << "unknown";
      break;
  }

  return os;
}

} // namespace core
