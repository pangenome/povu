#ifndef POVU_BOUNDED_Q_HPP
#define POVU_BOUNDED_Q_HPP

#include <condition_variable>
#include <deque>
#include <mutex>
#include <optional>
#include <utility>

namespace povu::bounded_queue {

// Multi-producer, multi-consumer bounded blocking queue (move-friendly).
template <class T>
class bounded_queue {
public:
  explicit bounded_queue(std::size_t capacity) : cap_(capacity ? capacity : 1) {}

  // Blocks when full; returns false if queue has been closed.
  bool push(T v) {
    std::unique_lock<std::mutex> lk(mx_);
    cv_not_full_.wait(lk, [&]{ return closed_ || q_.size() < cap_; });
    if (closed_) return false;
    q_.emplace_back(std::move(v));
    cv_not_empty_.notify_one();
    return true;
  }

  // Blocks when empty; returns std::nullopt if closed and drained.
  std::optional<T> pop() {
    std::unique_lock<std::mutex> lk(mx_);
    cv_not_empty_.wait(lk, [&]{ return closed_ || !q_.empty(); });
    if (q_.empty()) return std::nullopt;     // closed & drained
    T v = std::move(q_.front());
    q_.pop_front();
    cv_not_full_.notify_one();
    return v;
  }

  // Non-blocking helpers (optional)
  bool try_push(T v) {
    std::lock_guard<std::mutex> lk(mx_);
    if (closed_ || q_.size() >= cap_) return false;
    q_.emplace_back(std::move(v));
    cv_not_empty_.notify_one();
    return true;
  }

  std::optional<T> try_pop() {
    std::lock_guard<std::mutex> lk(mx_);
    if (q_.empty()) return std::nullopt;
    T v = std::move(q_.front());
    q_.pop_front();
    cv_not_full_.notify_one();
    return v;
  }

  void close() {
    std::lock_guard<std::mutex> lk(mx_);
    closed_ = true;
    cv_not_empty_.notify_all();
    cv_not_full_.notify_all();
  }

  bool closed() const {
    std::lock_guard<std::mutex> lk(mx_);
    return closed_;
  }

private:
  std::size_t cap_;
  std::deque<T> q_;
  mutable std::mutex mx_;
  std::condition_variable cv_not_empty_, cv_not_full_;
  bool closed_{false};
};

}


// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pbq = povu::bounded_queue;

#endif // POVU_BOUNDED_Q_HPP
