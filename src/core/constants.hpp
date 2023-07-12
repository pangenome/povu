#include <string>
#include <limits>



namespace core::constants {

  // colors
  const std::string gray{"gray"};
  const std::string black{"black"};
  const std::string red{"red"};

  // numeric
  const std::size_t SIZE_T_MAX = std::numeric_limits<size_t>::max();
  const std::size_t UNDEFINED_SIZE_T = std::numeric_limits<size_t>::max();
}
