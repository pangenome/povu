#include "./core.hpp"

namespace core {
std::ostream& operator<<(std::ostream& os, const task_t& t)  {
    switch (t) {
    case task_t::call:
    os << "call";
    break;
    default:
    os << "unknown";
    break;
    }
    return os;
}

} // namespace core
