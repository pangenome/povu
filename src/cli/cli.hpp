#include "./app.hpp"

namespace cli {

#define FILE_ERROR(name)                                                       \
  {                                                                            \
    std::string e = "Error, Failed to open the file " + name; \
    throw std::invalid_argument(e);                                              \
  }

// app version string constant
const std::string VERSION = "0.0.0-alpha";


int cli(int argc, char **argv, core::config& app_config);
} // namespace cli
