#ifndef POVU_COMPAT_HPP
#define POVU_COMPAT_HPP


#if __cplusplus >= 202002L && defined(__has_include)
  #if __has_include(<format>)
    #define HAS_STD_FORMAT 1
    #include <format>
  #endif
#endif

#ifndef HAS_STD_FORMAT
  #define HAS_STD_FORMAT 0
  #include <fmt/format.h>
#endif

namespace povu::compat {
#if HAS_STD_FORMAT
  using std::format;
#else
  using fmt::format;
#endif
} // namespace povu::compat

// add namespace alias for povu::compat
// below tells clang-tidy to skip that specific check for the next line.
// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pv_cmp = povu::compat;

#endif // POVU_COMPAT_HPP
