#include "./thread.hpp"

namespace povu::thread {

// Ceil division for non-negative integers
constexpr std::size_t ceil_div(std::size_t a, std::size_t b) {
  return b == 0 ? 0 : (a + b - 1) / b;
}

pt::op_t<u_int32_t> compute_thread_allocation(std::size_t requested_threads,
                                              std::size_t item_count) {
  // hw threads may be 0 on some platforms; clamp to at least 1
  std::size_t hw =
      std::max<std::size_t>(1, std::thread::hardware_concurrency());

  // If no items, run single-thread with chunk 0
  if (item_count == 0) {
    return {static_cast<u_int32_t>(1), static_cast<u_int32_t>(0)};
  }

  // Choose T: at least 1, at most hw, at most item_count
  std::size_t T =
    std::max<std::size_t>(1,std::min({requested_threads ? requested_threads : hw, hw, item_count}));

  // Ceil-divide so the last chunk isn't starved
  std::size_t chunk = ceil_div(item_count, T);

  return {static_cast<u_int32_t>(T), static_cast<u_int32_t>(chunk)};
}
} // namespace povu::thread
