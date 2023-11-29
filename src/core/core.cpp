#include "./core.hpp"

namespace core {

  std::ostream& operator<<(std::ostream& os, const color& c) {
		switch (c) {
		case color::gray:
			os << "gray";
			break;
		case color::black:
			os << "black";
			break;
		default:
			os << "unknown";
			break;
		}
		return os;
  }
  
  
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
