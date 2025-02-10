#include "./app.hpp"

namespace core {
std::ostream& operator<<(std::ostream& os, const task_e& t)  {
  switch (t) {
    case task_e::call:
      os << "call";
      break;
    case task_e::deconstruct:
      os << "deconstruct";
      break;
    case task_e::info:
      os << "info";
      break;
    default:
      os << "unknown";
      break;
  }

  return os;
}

} // namespace core
