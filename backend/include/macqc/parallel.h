// Thread helpers.
//
// run_parallel(threads, body) runs body(thread_index) for thread_index = 0 … threads − 1 at the
// same time and rethrows the first exception after all have finished. 23.0: the bodies run on a
// persistent pool of worker threads (created once, as many as the processor has hardware threads)
// instead of threads started and joined for every call, so that short parallel regions (a vector
// operation of an iterative solver, a level of a triangular solve) cost microseconds instead of a
// thread start each. All `threads` bodies run concurrently, as before (a body may wait for the
// others, e.g. at a SpinBarrier). When the pool is busy (two calculations at the same time, in
// two windows of the app) or a body calls run_parallel itself, the call starts its own threads,
// as 22.0 did.
//
// parallel_chunks(threads, chunks, f) runs f(c) for c = 0 … chunks − 1 on up to `threads`
// threads (static contiguous ranges of chunks). Reductions that sum one partial result per chunk
// in chunk order give the same bits for every thread count.
#pragma once

#include <algorithm>
#include <atomic>
#include <condition_variable>
#include <cstddef>
#include <exception>
#include <functional>
#include <mutex>
#include <thread>
#include <vector>

#include "macqc/job.h"

namespace macqc {

namespace detail {

// The pool of worker threads behind run_parallel. One parallel region at a time (`busy_`); a
// region gets the caller (index 0) and threads − 1 workers.
class ThreadPool {
public:
    static ThreadPool& instance() {
        static ThreadPool pool;
        return pool;
    }

    std::size_t workers() const { return threads_.size(); }

    // True if the region ran on the pool; false if the pool is busy or too small (the caller then
    // runs it another way).
    bool try_run(int threads, const std::function<void(int)>& body) {
        if (threads - 1 > static_cast<int>(threads_.size()) || in_worker()) return false;
        bool expected = false;
        if (!busy_.compare_exchange_strong(expected, true, std::memory_order_acquire)) return false;
        std::exception_ptr error;
        std::mutex error_mutex;
        Job* job = current_job();
        std::function<void(int)> wrapped = [&](int t) {
            JobScope scope(job);  // workers see the caller's job (cancellation, progress)
            try {
                body(t);
            } catch (...) {
                std::lock_guard<std::mutex> lock(error_mutex);
                if (!error) error = std::current_exception();
            }
        };
        {
            std::lock_guard<std::mutex> lock(mutex_);
            task_ = &wrapped;
            participants_ = threads - 1;
            remaining_.store(threads - 1, std::memory_order_relaxed);
            ++generation_;
        }
        wake_.notify_all();
        wrapped(0);
        // the workers finish their bodies; short regions are waited for by spinning
        for (int spin = 0; remaining_.load(std::memory_order_acquire) != 0; ++spin) {
            if (spin < 4096) continue;
            std::unique_lock<std::mutex> lock(mutex_);
            done_.wait(lock, [&] { return remaining_.load(std::memory_order_acquire) == 0; });
            break;
        }
        {
            std::lock_guard<std::mutex> lock(mutex_);
            task_ = nullptr;
        }
        busy_.store(false, std::memory_order_release);
        if (error) std::rethrow_exception(error);
        return true;
    }

    ~ThreadPool() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            stop_ = true;
            ++generation_;
        }
        wake_.notify_all();
        for (auto& t : threads_) t.join();
    }

    ThreadPool(const ThreadPool&) = delete;
    ThreadPool& operator=(const ThreadPool&) = delete;

private:
    ThreadPool() {
        const unsigned hw = std::thread::hardware_concurrency();
        const int n = hw > 1 ? static_cast<int>(hw) - 1 : 0;
        for (int i = 0; i < n; ++i) threads_.emplace_back([this, i] { loop(i + 1); });
    }

    static bool& in_worker() {
        static thread_local bool flag = false;
        return flag;
    }

    void loop(int index) {
        in_worker() = true;
        unsigned long long seen = 0;
        for (;;) {
            const std::function<void(int)>* task = nullptr;
            {
                std::unique_lock<std::mutex> lock(mutex_);
                wake_.wait(lock, [&] { return stop_ || generation_ != seen; });
                if (stop_) return;
                seen = generation_;
                if (index > participants_) continue;  // not needed for this region
                task = task_;
            }
            (*task)(index);
            if (remaining_.fetch_sub(1, std::memory_order_acq_rel) == 1) {
                std::lock_guard<std::mutex> lock(mutex_);
                done_.notify_all();
            }
        }
    }

    std::vector<std::thread> threads_;
    std::mutex mutex_;
    std::condition_variable wake_, done_;
    const std::function<void(int)>* task_ = nullptr;
    int participants_ = 0;
    unsigned long long generation_ = 0;
    bool stop_ = false;
    std::atomic<int> remaining_{0};
    std::atomic<bool> busy_{false};
};

}  // namespace detail

template <class Body>
void run_parallel(int threads, Body body) {
    if (threads <= 1) {
        body(0);
        return;
    }
    const std::function<void(int)> f = [&body](int t) { body(t); };
    if (detail::ThreadPool::instance().try_run(threads, f)) return;
    // the pool is busy or too small: threads of this call (22.0)
    std::vector<std::thread> pool;
    std::exception_ptr error;
    std::mutex mutex;
    Job* job = current_job();  // workers see the caller's job (cancellation, progress)
    for (int t = 0; t < threads; ++t)
        pool.emplace_back([&, t, job] {
            JobScope scope(job);
            try {
                body(t);
            } catch (...) {
                std::lock_guard<std::mutex> lock(mutex);
                if (!error) error = std::current_exception();
            }
        });
    for (auto& th : pool) th.join();
    if (error) std::rethrow_exception(error);
}

// f(c) for c in [0, chunks): thread t of `threads` takes the contiguous range of chunks
// [chunks·t/threads, chunks·(t+1)/threads).
template <class F>
void parallel_chunks(int threads, std::size_t chunks, F f) {
    const int t = static_cast<int>(std::min<std::size_t>(static_cast<std::size_t>(std::max(1, threads)), chunks));
    if (t <= 1) {
        for (std::size_t c = 0; c < chunks; ++c) f(c);
        return;
    }
    run_parallel(t, [&](int k) {
        const std::size_t a = chunks * static_cast<std::size_t>(k) / static_cast<std::size_t>(t);
        const std::size_t b = chunks * static_cast<std::size_t>(k + 1) / static_cast<std::size_t>(t);
        for (std::size_t c = a; c < b; ++c) f(c);
    });
}

// A barrier for the `count` concurrent bodies of one run_parallel region (sense-reversing; spins
// with the processor's pause hint, after 65536 rounds it yields the processor).
class SpinBarrier {
public:
    explicit SpinBarrier(int count) : count_(count), waiting_(0), phase_(0) {}

    void wait() {
        const unsigned phase = phase_.load(std::memory_order_acquire);
        if (waiting_.fetch_add(1, std::memory_order_acq_rel) + 1 == count_) {
            waiting_.store(0, std::memory_order_relaxed);
            phase_.store(phase + 1, std::memory_order_release);
            return;
        }
        for (int spin = 0; phase_.load(std::memory_order_acquire) == phase; ++spin) {
            if (spin < (1 << 16)) {
                cpu_relax();
            } else {
                std::this_thread::yield();
            }
        }
    }

private:
    static void cpu_relax() {
#if defined(__x86_64__) || defined(__i386__)
        __builtin_ia32_pause();
#elif defined(__aarch64__)
        __asm__ __volatile__("yield");
#endif
    }

    const int count_;
    std::atomic<int> waiting_;
    std::atomic<unsigned> phase_;
};

}  // namespace macqc
