#ifndef POVU_LOG_HPP
#define POVU_LOG_HPP

#include <fmt/format.h>
#include <fmt/color.h>


#define FN() pv_cmp::format("[{}::{}]", MODULE, __func__)

#define DEBUG_PRINT(format, ...)                                               \
  fprintf(stderr, "%s:%d: " format, __FILE__, __LINE__, ##__VA_ARGS__)

/**
 * Generic logging macro: prints a label (e.g. "ERR") with color and formatted
 * output.
 */
#define LOG(label, color, fmt_str, ...)                                        \
  do {                                                                         \
    fmt::print(stderr, fmt::fg(color) | fmt::emphasis::bold, "{} {} ", label,  \
               FN());                                                          \
    fmt::print(stderr, fmt_str, ##__VA_ARGS__);                                \
    fmt::print(stderr, "\n");                                                  \
  } while (false)

#define ERR(fmt_str, ...)                                                      \
  LOG("ERR", fmt::color::crimson, fmt_str, ##__VA_ARGS__)

#define WARN(fmt_str, ...)                                                     \
  LOG("WARN", fmt::color::yellow, fmt_str, ##__VA_ARGS__)


#endif // POVU_LOG_HPP
