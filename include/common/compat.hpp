#ifndef POVU_COMPAT_HPP
#define POVU_COMPAT_HPP


// Required for erase_if fallback
#if __cplusplus < 202002L
  #include <algorithm>
  #include <iterator>
#endif

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

// erase_if
// erase_if is C++20
// remove elements that do not meet condition iv and v
  // std::erase_if(src_lca_vec, not_cond_ii_iii);

  // remove_if is C++11
  //   std::remove_if(...) moves the elements that don't match the predicate to
  //   the end. .erase(...) trims those elements off the container.

#if __cplusplus >= 202002L
  using std::erase_if;
#else
  template <typename Container, typename Predicate>
  void erase_if(Container &c, Predicate pred) {
    c.erase(std::remove_if(c.begin(), c.end(), pred), c.end());
  }
#endif

// backport contains for C++11 and C++14
template <typename Container, typename Key>
bool contains(const Container &c, const Key &key) {
  return c.find(key) != c.end();
}

} // namespace povu::compat

// add namespace alias for povu::compat
// below tells clang-tidy to skip that specific check for the next line.
// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pv_cmp = povu::compat;

#endif // POVU_COMPAT_HPP
