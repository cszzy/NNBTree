#include <mutex>
#include <atomic>
#include <vector>
#include <thread>
#include <iostream>

class __attribute__((aligned(64))) Spinlock {
  public:
  Spinlock() = default;
  Spinlock(const Spinlock&) = delete;
  Spinlock& operator= (const Spinlock&) = delete;
  void lock() {
      // bool lock = false;
      while (flag.test_and_set(std::memory_order_acquire)) {
          // if (log_lock) {
          //     if (lock == false) {
          //         lock = true;
          //         lock_times++;
          //     }
          // }
      }
  }

  void unlock() { flag.clear(std::memory_order_release); }

  private:
  std::atomic_flag flag = ATOMIC_FLAG_INIT;
  char padding[63];
};

int main() {
    int thread_num = 4;
    const int lock_nums = 10;
    std::vector<Spinlock*> lock_space; 
    for (int i = 0; i < lock_nums; i++) {
        lock_space.emplace_back(new Spinlock());
    }
    std::vector<std::thread> threads;
    for(int i = 0; i < thread_num; i ++) {
        threads.emplace_back([&](){
            int i = rand() % lock_nums;
            lock_space[i]->lock();
            for (int j = 0; j < 10000000; j++) {
                volatile int k = 10;
                k *= j;
            }
            lock_space[i]->unlock();
        });
    }

    for (auto &t : threads) {
        t.join();
    }

    return 0;
}