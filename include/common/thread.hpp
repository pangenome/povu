#ifndef POVU_THREAD_HPP
#define POVU_THREAD_HPP

#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <exception>
#include <functional>
#include <mutex>
#include <queue>
#include <thread>
#include <vector>
#include <future>
#include <type_traits>
#include <utility>

#include "./compat.hpp"
#include "./core.hpp"


namespace povu::thread {

class thread_pool {
public:
    explicit thread_pool(std::size_t n) : stop_(false) {
        if (n == 0) n = 1;
        workers_.reserve(n);
        for (std::size_t i = 0; i < n; ++i)
            workers_.emplace_back([this]{ worker(); });
    }
    ~thread_pool() {
        { std::lock_guard<std::mutex> lk(mx_); stop_ = true; }
        cv_.notify_all();
        for (auto& t : workers_) t.join();
    }

    std::size_t size() const noexcept { return workers_.size(); }

    // fire-and-forget tasks
    template <class F>
    void enqueue(F&& f) {
        {
            std::lock_guard<std::mutex> lk(mx_);
            jobs_.emplace(std::forward<F>(f));
        }
        cv_.notify_one();
    }


    template <class F>
    auto submit(F &&f) -> std::future<std::invoke_result_t<F &>>;

  private:
    void worker() {
        for (;;) {
            std::function<void()> job;
            {
                std::unique_lock<std::mutex> lk(mx_);
                cv_.wait(lk, [&]{ return stop_ || !jobs_.empty(); });
                if (stop_ && jobs_.empty()) return;
                job = std::move(jobs_.front());
                jobs_.pop();
            }
            job(); // run outside lock
        }
    }

    std::vector<std::thread> workers_;
    std::queue<std::function<void()>> jobs_;
    std::mutex mx_;
    std::condition_variable cv_;
    bool stop_;
};

// A small "task group" to wait for a set of enqueued tasks (no futures)
class task_group {
public:
    explicit task_group(thread_pool& p) : pool_(p), count_(0) {}

    template <class F>
    void run(F&& f) {
        count_.fetch_add(1, std::memory_order_relaxed);
        pool_.enqueue([this, func = std::forward<F>(f)]() mutable {
            try { func(); }
            catch (...) {
                std::lock_guard<std::mutex> lk(ex_mx_);
                errors_.push_back(std::current_exception());
            }
            if (count_.fetch_sub(1, std::memory_order_acq_rel) == 1) {
                std::lock_guard<std::mutex> lk(wait_mx_);
                wait_cv_.notify_all();
            }
        });
    }

    void wait() {
        std::unique_lock<std::mutex> lk(wait_mx_);
        wait_cv_.wait(lk, [&]{ return count_.load(std::memory_order_acquire) == 0; });
        if (!errors_.empty()) std::rethrow_exception(errors_.front()); // or handle all
    }

private:
    thread_pool& pool_;
    std::atomic<std::size_t> count_;
    std::mutex wait_mx_;
    std::condition_variable wait_cv_;

    std::mutex ex_mx_;
    std::vector<std::exception_ptr> errors_;
};

// Template definitions must be visible to callers.
template <class F>
auto thread_pool::submit(F &&f) -> std::future<std::invoke_result_t<F &>> {
  using R = std::invoke_result_t<F &>;

  // Make the queued callable copyable for std::function
  auto task = std::make_shared<std::packaged_task<R()>>(std::forward<F>(f));
  std::future<R> fut = task->get_future();

  {
    std::lock_guard<std::mutex> lk(mx_);
    jobs_.emplace([task] { (*task)(); });
  }
  cv_.notify_one();
  return fut;
}


/**
 * Divide the number of components into chunks for each thread
 * @param tc: (thread_count) number of threads to use
 * @param ic: (item_count) number of items to process
 */
// std::pair<uint32_t, uint32_t>
pt::op_t<u_int32_t> compute_thread_allocation(std::size_t requested_threads,
                                              std::size_t item_count);

} // namespace povu::thread

#endif
