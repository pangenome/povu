#ifndef UTILS_HPP
#define UTILS_HPP

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <cstddef>
#include <ctime>
#include <functional>
#include <future>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <mutex>
#include <queue>
#include <sstream>
#include <string>
#include <sys/types.h>
#include <thread>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include "./compat.hpp"
#include "./types/core.hpp"


namespace povu::utils {

std::string reverse_complement(const std::string &sequence);

/**
 * @brief returns true if a string is made up of only digits.
 */
bool is_numeric_string(const std::string &s);


template <typename Container> std::string concat_with(const Container &v, char delim) {
  std::ostringstream oss;
  auto it = v.begin();
  if (it != v.end()) {
    oss << *it;
    ++it;
    for (; it != v.end(); ++it) {
      oss << delim << *it;
    }
  }
  return oss.str();
}

// TODO rename to print_with_delim or print_with
template <typename T> void print_with_comma(std::ostream& os, const T& v, char delim) {
  if (v.empty()) { return; }

  for (auto it {v.begin()}; it != v.end(); ++it) {
    os << *it;
    if (std::next(it) != v.end()) { os << delim << " "; }
  }
}

/**
  * @brief
 */
void report_time(std::ostream& os, std::string fn_name, std::string action, std::chrono::duration<double> period);


/**
 * @brief Returns the current date in the format YYYYMMDD
 */
std::string today();

/**
 * @brief
 *
 * @param v: the vector whose value is to be erased passed by copy to avoid mutating the original
 * @param idx: the index to be erased
 * @return a vector with the value at the given index erased
 */
std::vector<std::string> immutable_erase(std::vector<std::string> v, std::size_t idx);

// TODO : move to povu::types
template <typename Key, typename Value>
class TwoWayMap {
    std::unordered_map<Key, Value> keyToValueMap;
    std::unordered_map<Value, Key> valueToKeyMap;

public:
    // Insert a key-value pair
    void insert(const Key& key, const Value& value) {
        keyToValueMap[key] = value;
        valueToKeyMap[value] = key;
    }

    // Lookup value by key
    Value get_value(const Key& key) const {
        auto it = keyToValueMap.find(key);
        if (it != keyToValueMap.end()) {
            return it->second;
        }
        // Return a default-constructed Value or throw an exception, depending on your use case
        return Value{};
    }

    // Lookup key by value
    Key get_key(const Value& value) const {
        auto it = valueToKeyMap.find(value);
        if (it != valueToKeyMap.end()) {
            return it->second;
        }
        // Return a default-constructed Key or throw an exception, depending on your use case
        return Key{};
    }

  // get all keys
  std::vector<Key> get_keys() const {
    std::vector<Key> keys;
    for (const auto& [key, _]: keyToValueMap) {
      keys.push_back(key);
    }
    return keys;
  }

  // find key
  bool has_key(const Key& key) const {
    return keyToValueMap.find(key) != keyToValueMap.end();
  }
};

template <typename T> void push_front(std::vector<T>& v, const T& elem) {
  v.insert(v.begin(), elem);
}

/**
 * @brief given a bi-Edged index return its index in the bi-Directed graph
 *
 * @param
 * @param
 * @return
 */
std::size_t to_bidirected_idx(std::size_t be_idx, bool has_dummy=true);

/**
 * @brief given a bi-Directed index return the indexes of the left and right vertices in the bi-Edged graph
 *
 * @param
 * @param
 * @return
 */
std::pair<std::size_t, std::size_t> frm_bidirected_idx(std::size_t bd_idx, bool has_dummy=true);

/**
 * @brief split a string into tokens using a delimiter
 */
void split(const std::string &line, char sep, std::vector<std::string>* tokens);

inline thread_local bool tp_in_worker = false;

class ThreadPool {
public:
  explicit ThreadPool(std::size_t threads = std::thread::hardware_concurrency())
      : stop_(false) {
    if (threads == 0) threads = 1;
    workers_.reserve(threads);
    for (std::size_t i = 0; i < threads; ++i) {
      workers_.emplace_back([this] {
        tp_in_worker = true;
        for (;;) {
          Task task;
          { // wait for work or shutdown
            std::unique_lock<std::mutex> lk(mx_);
            cv_.wait(lk, [this] { return stop_ || !q_.empty(); });
            if (stop_ && q_.empty()) return;
            task = std::move(q_.front());
            q_.pop();
          }
          task(); // run outside the lock
        }
      });
    }
  }

  // Non-copyable, movable if you like (omitted here for brevity)
  ThreadPool(const ThreadPool&) = delete;
  ThreadPool& operator=(const ThreadPool&) = delete;

  ~ThreadPool() {
    { std::lock_guard<std::mutex> lk(mx_); stop_ = true; }
    cv_.notify_all();
    for (auto& t : workers_) t.join();
  }

  template <class F, class... Args>
  auto enqueue(F&& f, Args&&... args)
      -> std::future<std::invoke_result_t<F, Args...>> {
    using R = std::invoke_result_t<F, Args...>;
    auto task_ptr = std::make_shared<std::packaged_task<R()>>(
        std::bind(std::forward<F>(f), std::forward<Args>(args)...));
    std::future<R> fut = task_ptr->get_future();
    {
      std::lock_guard<std::mutex> lk(mx_);
      if (stop_) throw std::runtime_error("ThreadPool stopped");
      q_.emplace([task_ptr] { (*task_ptr)(); });
    }
    cv_.notify_one();
    return fut;
  }

  std::size_t size() const noexcept { return workers_.size(); }

private:
  using Task = std::function<void()>;
  std::vector<std::thread> workers_;
  std::queue<Task> q_;
  std::mutex mx_;
  std::condition_variable cv_;
  bool stop_;
};


// Assisted nested parallel-for: caller participates; helpers = T or T-1 if already inside pool
template <class F>
void parallel_for_assisted(ThreadPool& pool, std::size_t N, F&& body, std::size_t block = 0) {
  if (N == 0) return;

  const std::size_t T = std::max<std::size_t>(1, pool.size());
  if (block == 0) block = std::max<std::size_t>(1, N / (T * 16));

  std::atomic<std::size_t> next{0};

  // If already running inside a pool worker, keep one worker "free": spawn T-1 helpers.
  const std::size_t helpers = tp_in_worker ? (T > 1 ? T - 1 : 0) : T;

  std::vector<std::future<void>> fs;
  fs.reserve(helpers);
  for (std::size_t h = 0; h < helpers; ++h) {
    fs.emplace_back(pool.enqueue([&, block] {
      for (;;) {
        const std::size_t start = next.fetch_add(block, std::memory_order_relaxed);
        if (start >= N) break;
        const std::size_t end = std::min(start + block, N);
        for (std::size_t i = start; i < end; ++i) body(i);
      }
    }));
  }

  // Caller helps execute the same loop (crucial for nested calls).
  for (;;) {
    const std::size_t start = next.fetch_add(block, std::memory_order_relaxed);
    if (start >= N) break;
    const std::size_t end = std::min(start + block, N);
    for (std::size_t i = start; i < end; ++i) body(i);
  }

  for (auto& f : fs) f.get(); // rethrow exceptions
}

/**
 * Divide the number of components into chunks for each thread
 * @param tc: (thread_count) number of threads to use
 * @param ic: (item_count) number of items to process
 */
// std::pair<uint32_t, uint32_t>
pt::op_t<u_int32_t> compute_thread_allocation(std::size_t requested_threads,
                                              std::size_t item_count);

} // namespace povu::utils


// NOLINTNEXTLINE(misc-unused-alias-decls)
namespace pu = povu::utils;
#endif
