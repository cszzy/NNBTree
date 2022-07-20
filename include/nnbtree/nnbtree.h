#pragma once

#include <cassert>
#include <climits>
#include <fstream>
#include <future>
#include <iostream>
#include <math.h>
#include <mutex>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#include <vector>
#include <unordered_set>

#include "nvm_alloc.h"
#include "common.h"
#include "numa_config.h"
#include "tree_log.h"

#define USE_SPINLOCK // 使用spinlock还是mutex

#define CACHE_SUBTREE // 是否开启子树缓存

#define MIGRATE_SUBTREE // 是否开启子树迁移

#define BG_THREAD_NUMS 1

#define PAGESIZE 256

#define CACHE_LINE_SIZE 64

#define TOPK_SUBTREE_NUM 1000

// 指示lookup移动方向
#define IS_FORWARD(c) (c % 2 == 0)

#define MAX_SUBTREE_HEIGHT 3 // 子树最大高度, 应尽量地设置小一些提高性能

using entry_key_t = uint64_t;

namespace nnbtree {

class Spinlock {
  public:
  Spinlock() : flag(ATOMIC_FLAG_INIT) {};
  Spinlock(const Spinlock&) = delete;
  Spinlock& operator= (const Spinlock&) = delete;
  ~Spinlock() {std::cout << "spinlock " << this << "destructed" << std::endl << std::flush; }
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
  std::atomic_flag flag{ATOMIC_FLAG_INIT};
  // char padding[63];
};

// const size_t NVM_ValueSize = 256;
static inline void alloc_memalign(void **ret, size_t alignment, size_t size) {
#ifdef USE_MEM
  int res = posix_memalign(ret, alignment, size);
  p_assert(res == 0, "mem alloc fail");
#else
  // *ret =  NVM::data_alloc->alloc(size);
  *ret = nnbtree::index_pmem_alloc(size);
  p_assert(*ret, "mem/pmem alloc fail");
#endif
}

class Page;

enum SubTreeStatus : uint32_t {
  IN_NVM,
  IN_DRAM,
  NEED_MOVE_TO_NVM,
  NEED_MOVE_TO_DRAM, // 迁移到
  IS_DELETED, // 防止迁移已经完全删除的子树
  NEED_MOVE_TO_REMOTE // 迁移到远程
};

class SubTree {
private:
  int32_t height_; // 4B
  // int32_t numa_id; // 标记子树在哪个numa?但是子树下的page仍然处于不同的numa节点=_=
  SubTreeStatus subtree_status_; // 标识subtree状态
  // char *root; // 整个B+树的根
  char *sub_root_; // 8B 正在使用的subtree的根
  char *nvm_root_; // 8B 当子树移到内存时，用此记录nvm的root
  SubTree *left_sibling_subtree_; //左子树
  SubTree *right_sibling_subtree_; //右子树
  uint64_t read_times_[numa_node_num]; // 各numa节点写操作次数
  uint64_t write_times_[numa_node_num]; // 各numa节点读操作次数
  TreeLog *treelog_; // 写日志
  uint64_t hotness_[numa_node_num]; // 记录的是上一周期的热度
  uint64_t tmp_hotness; // 只使用缓存版本，记录当前热度
  
  // bool is_dirty_; // 脏页写回NVM
#ifndef USE_SPINLOCK
  std::mutex subtree_lock; // 每个子树一个锁
#else
  Spinlock subtree_lock;
#endif

public:
  SubTree(SubTreeStatus status = SubTreeStatus::IN_NVM);
  void setNewRoot(char *); // 设置sub_root
  void getNumberOfNodes();
  void btree_insert(entry_key_t, char *);
  void btree_insert_internal(char *, entry_key_t, char *, uint32_t);
  void btree_delete(entry_key_t);
  void btree_delete_internal(entry_key_t, char *, uint32_t, entry_key_t *,
                             bool *, Page **);
  char *btree_search(entry_key_t);
  void btree_search_range(entry_key_t, entry_key_t, unsigned long *);
  void printAll();

  // zzy add
  SubTree(Page *root, SubTreeStatus status = SubTreeStatus::IN_NVM);
  Page *btree_search_leaf(entry_key_t);
  void btree_search_range(entry_key_t, entry_key_t, std::vector<std::pair<uint64_t, uint64_t>> &result, int &size); 
  void btree_search_range(entry_key_t, entry_key_t, void **values, int &size); 
  void PrintInfo();
  void CalculateSapce(uint64_t &space);
  Page *getRoot();

  void cal_hotness();
  uint64_t get_hotness() const;
  void setSubTreeStatus(SubTreeStatus status);
  SubTreeStatus getSubTreeStatus() const;
  void move_to_nvm(); // 淘汰
  void move_to_dram(); // 缓存

  void lock_subtree();
  void unlock_subtree();

  friend class Page;
};

class IndexTree {
private:
  int height_;
  char *root; // 整个B+树的根

public:
  IndexTree(Page*, uint32_t);
  void setNewRoot(char *); // 设置root
  void getNumberOfNodes();
  void btree_insert(entry_key_t, char *);
  void btree_insert_internal(char *, entry_key_t, char *, uint32_t);
  void btree_delete(entry_key_t);
  void btree_delete_internal(entry_key_t, char *, uint32_t, entry_key_t *,
                             bool *, Page **);
  char *btree_search(entry_key_t);
  void btree_search_range(entry_key_t, entry_key_t, unsigned long *);
  void printAll();

  // zzy add
  // Page *btree_search_leaf(entry_key_t);
  void btree_search_range(entry_key_t, entry_key_t, std::vector<std::pair<uint64_t, uint64_t>> &result, int &size); 
  void btree_search_range(entry_key_t, entry_key_t, void **values, int &size); 
  void PrintInfo();
  void CalculateSapce(uint64_t &space);
  Page *getRoot();

  friend class Page;
};

class NNBTree;
NNBTree *index_tree_root;

class Statistics {
    public:
        Statistics() {}

        ~Statistics() {}

        void insert_subtree(SubTree * subtree) {
            all_subtree_lock.lock();
            all_subtree.push_back(subtree);
            all_subtree_lock.unlock();
        }

        // 选择前k个热点子树，目前只考虑缓存
        void select_topk(int k) {  
            // 把子树复制过来 O(N/k)
            all_subtree_lock.lock();
            staticstic_subtree.insert(staticstic_subtree.end(), 
                    all_subtree.begin() + staticstic_subtree.size(), all_subtree.end());
            all_subtree_lock.unlock();

            if (unlikely(staticstic_subtree.size() < k)) { // 子树总数小于缓存数
                return;
            }

            // 计算所有子树热度 O(N)
            // 这一步现在放到topk中进行计算，减少一次扫描
            // for (auto subtree : staticstic_subtree) {
            //     subtree->cal_hotness();
            // }

            // O(N)得到topk，同时将热点nvm子树设置为待缓存，将非热点nvm子树设置为待淘汰
            topk(staticstic_subtree, k);
            std::cout << "select topk" << std::endl;
        }

    private:
        void topk(std::vector<SubTree *> &arr, int k) {
            if (unlikely(arr.size() < k)) {
                // subtree总数小于k，不进行缓存
                return;
            }

            srand((unsigned)time(NULL));
            int left = 0;
            int right = arr.size() - 1;
            int target = arr.size() - k;
            while (true) {
                int pivotIndex = partition(arr, left, right);
                if (pivotIndex == target) {
                    // 设置前k个nvm子树的状态为NEED_MOVE_TO_DRAM
                    for (int i = arr.size() - 1; i >= target; i--) {
                        topk_subtree.erase(arr[i]);
                        if (arr[i]->getSubTreeStatus() == SubTreeStatus::IN_NVM) {
                            arr[i]->setSubTreeStatus(SubTreeStatus::NEED_MOVE_TO_DRAM);
                        }
                        //  else if (arr[i]->getSubTreeStatus() == SubTreeStatus::NEED_MOVE_TO_NVM) {
                        //     arr[i]->setSubTreeStatus(SubTreeStatus::IN_DRAM);
                        //     // std::cout << "dram subtree hasn't been move to nvm" << std::endl;
                        // }
                    }

                    for (auto iter = topk_subtree.begin(); iter != topk_subtree.end(); iter++) {
                        if ((*iter)->getSubTreeStatus() == SubTreeStatus::IN_DRAM) {
                                (*iter)->setSubTreeStatus(SubTreeStatus::NEED_MOVE_TO_NVM);
                            }
                    }

                    // 重置topk
                    topk_subtree.clear();
                    for (int i = arr.size() - 1; i >= target; i--) {
                        topk_subtree.insert(arr[i]);
                    }

                    return; 
                } else if (pivotIndex < target) {
                    left = pivotIndex + 1; 
                } else {
                    // pivotIndex > target
                    right = pivotIndex - 1; 
                }
            }
        }

        int partition(std::vector<SubTree *> &arr, int low, int high) {
            bool need_cal_hotness = false;
            if (low == 0 && high == arr.size() - 1) {
                // 第一次partition，需要计算热度
                need_cal_hotness = true;
            }
            int pos = low + random() % (high - low + 1);
            // int pos = (low + high) >> 1;
            if (need_cal_hotness) {
                arr[pos]->cal_hotness();
            }
            uint64_t pivot = arr[pos]->get_hotness();
            std::swap(arr[low], arr[pos]);
            int left = low + 1;
            int right = high;
            if (need_cal_hotness) {
                arr[left]->cal_hotness();
                arr[right]->cal_hotness();
            }

            while (true) {
                while (left <= right && arr[left]->get_hotness() <= pivot) {
                    left++;
                    if (left <= right && need_cal_hotness) {
                        arr[left]->cal_hotness();
                    }
                }

                while (left <= right && arr[right]->get_hotness() >= pivot) {
                    right--;
                    if (left <= right && need_cal_hotness) {
                        arr[right]->cal_hotness();
                    }
                }

                if (left < right) {
                    std::swap(arr[left], arr[right]);
                    left++;
                    right--;
                } else {
                    break;
                }
            }

            std::swap(arr[low], arr[right]);

            return right;
        }

        std::vector<SubTree *> all_subtree; // 记录所有subtree指针
        std::mutex all_subtree_lock; // all_subtree锁，保护前台线程写和后线程读
        std::vector<SubTree *> staticstic_subtree; // 使用topk进行统计
        std::unordered_set<SubTree *> topk_subtree; // 记录topk子树
};

Statistics *statis_;

std::vector<std::thread> bg_thread; // 后台线程：用来统计和压缩

void bgthread_func(int bg_thread_id) {
    std::this_thread::sleep_for(std::chrono::seconds(3));
    while (true) {
        // statis_->select_topk(TOPK_SUBTREE_NUM);
        std::this_thread::sleep_for(std::chrono::seconds(15));
    }
}

void start_bgthread() {
    int bg_thread_num = 1;
    for (int i = 0; i < bg_thread_num; i++) {
        bg_thread.emplace_back(std::thread(bgthread_func, i));
    }
}

class NNBTree {
public:
    NNBTree() {
      has_indextree = false;
      tree_ = (char *)new SubTree();

      statis_ = new Statistics();
      statis_->insert_subtree((SubTree *)tree_);

      index_tree_root = this;
      std::cout << "index_tree_root: " << index_tree_root << std::endl;
    
      start_bgthread();
    }

    ~NNBTree() {
        delete tree_;
    }

    void set_indextree(const IndexTree * indextree_) {
        tree_ = (char *)indextree_;
        has_indextree = true;
    }

    bool has_hasindextree() {
        return has_indextree;
    }

    char *btree_search(entry_key_t key) {
        if (likely(has_indextree))
            return ((IndexTree *)tree_)->btree_search(key);
        return ((SubTree *)tree_)->btree_search(key);
    }

    void btree_insert(entry_key_t key, char *right) { // need to be string
        if (likely(has_indextree))
            return ((IndexTree *)tree_)->btree_insert(key, right);
        return ((SubTree *)tree_)->btree_insert(key, right);
    }

    void btree_delete(entry_key_t key) {
        if (likely(has_indextree))
            return ((IndexTree *)tree_)->btree_delete(key);
        return ((SubTree *)tree_)->btree_delete(key);
    }

    // Function to search keys from "min" to "max"
    void btree_search_range(entry_key_t min, entry_key_t max,
                                unsigned long *buf) {
        if (likely(has_indextree))
            return ((IndexTree *)tree_)->btree_search_range(min, max, buf);
        return ((SubTree *)tree_)->btree_search_range(min, max, buf);
    }

    void btree_search_range(entry_key_t min, entry_key_t max, 
        std::vector<std::pair<uint64_t, uint64_t>> &result, int &size) {
        if (likely(has_indextree))
            return ((IndexTree *)tree_)->btree_search_range(min, max, result, size);
        return ((SubTree *)tree_)->btree_search_range(min, max, result, size);
    }

    void btree_search_range(entry_key_t min, entry_key_t max, void **values, int &size) {
        if (likely(has_indextree))
            return ((IndexTree *)tree_)->btree_search_range(min, max, values, size);
        return ((SubTree *)tree_)->btree_search_range(min, max, values, size);
    }


    void PrintInfo() {
        if (likely(has_indextree))
            return ((IndexTree *)tree_)->PrintInfo();
        return ((SubTree *)tree_)->PrintInfo();
    }

    void printAll() {
        if (likely(has_indextree))
            return ((IndexTree *)tree_)->printAll();
        return ((SubTree *)tree_)->printAll();
    }

    void indextree_lock() {
      tree_lock.lock();
    }

    void indextree_unlock() {
      tree_lock.unlock();
    }

public:
    char *tree_;
    bool has_indextree;
    std::mutex tree_lock;
};

enum PageType : uint8_t {
  NVM_SUBTREE_PAGE,
  DRAM_CACHETREE_PAGE,
  DRAM_INDEXTREE_PAGE,
  INDEXTREE_LAST_LEVEL_PAGE, // 标识indextree最后一层page, 其内部存储的是subtree的指针
  UNKNOWN
};

class header {
private:
  Page *leftmost_ptr;     // 8 bytes // 非叶子节点最左侧子节点的指针
  Page *right_sibling_ptr;      // 8 bytes // 类似blink Tree，右兄弟结点指针
  // Page *left_sibling_ptr;  // 8 bytes // 类似blink Tree，右兄弟结点指针
#ifndef USE_SPINLOCK
  std::mutex *mtx;        // 8 bytes // 写独占
#else
  Spinlock *mtx; // subtree实际不使用这个字段
#endif
  uint32_t level;         // 4 bytes // 叶节点level为0, 向上累加
  uint8_t switch_counter; // 1 bytes // 指导读线程的扫描方向, 偶数代表为insert, 奇数代表delete
  uint8_t is_deleted;     // 1 bytes // ?实际没有节点合并, 用不到
  int8_t last_index;      // 1 bytes // 指示最后条目的位置
  uint8_t page_type;      // 1 bytes // 指示page类型, 0: nvm_subtree_page, 1: dram_cachetree_page, 2: dram_indextree_page
  // char padding[8]; // 填充
  

  friend class Page;
  friend class SubTree;
  friend class IndexTree;

public:
  header() : leftmost_ptr(nullptr), right_sibling_ptr(nullptr), mtx(nullptr), 
            level(0), switch_counter(0), is_deleted(false),
            last_index(-1), page_type(PageType::UNKNOWN) {
    // mtx = new std::mutex();
    // mtx = new Spinlock();
    // std::cout << " cons " << mtx << std::endl << std::flush;
  }

  ~header() { 
    if(mtx) {
      delete mtx; 
      mtx = NULL;
    } 
  }
};

class entry {
private:
  entry_key_t key; // 8 bytes
  char *ptr;       // 8 bytes

public:
  entry() {
    key = LONG_MAX;
    ptr = NULL;
  }

  friend class Page;
  friend class SubTree;
  friend class IndexTree;
};

// 计算page和cacheline可以存储的条目数目

const int cardinality = (PAGESIZE - sizeof(header)) / sizeof(entry);
const int count_in_line = CACHE_LINE_SIZE / sizeof(entry);

// B+树节点
class Page {
private:
  header hdr;                 // header in persistent memory, 32 bytes
  entry records[cardinality]; // slots in persistent memory, 16 bytes * n

public:
  friend class SubTree;
  friend class IndexTree;

  Page(PageType page_type_, uint32_t level) {
    p_assert(sizeof(Page) == PAGESIZE, "class Page size is not %d", PAGE_SIZE);
    hdr.level = level;
    records[0].ptr = NULL;
    hdr.page_type = (uint8_t)page_type_;
    if (page_type_ == PageType::DRAM_INDEXTREE_PAGE || page_type_ == PageType::INDEXTREE_LAST_LEVEL_PAGE) {
#ifndef USE_SPINLOCK
      hdr.mtx = new std::mutex();
#else
      hdr.mtx = new Spinlock();
#endif
    }
  }

  // this is called when tree grows
  Page(PageType page_type_, Page *left, entry_key_t key, Page *right, uint32_t level = 0) {
    p_assert(sizeof(Page) == PAGESIZE, "class Page size is not %d", PAGE_SIZE);
    hdr.leftmost_ptr = left;
    hdr.level = level;
    records[0].key = key;
    records[0].ptr = (char *)right;
    records[1].ptr = NULL;

    hdr.last_index = 0;
    hdr.page_type = (uint8_t)page_type_;
    if (page_type_ == PageType::DRAM_INDEXTREE_PAGE || page_type_ == PageType::INDEXTREE_LAST_LEVEL_PAGE) {
#ifndef USE_SPINLOCK
      hdr.mtx = new std::mutex();
#else
      hdr.mtx = new Spinlock();
#endif
    }

    if (page_is_inpmem()) {
      clflush((char *)this, sizeof(Page));
    }
  }

  // this is called when add new subtree
  Page(PageType page_type_, SubTree *left, entry_key_t key, SubTree *right, uint32_t level = 0) {
    p_assert(sizeof(Page) == PAGESIZE, "class Page size is not %d", PAGE_SIZE);
    hdr.leftmost_ptr = (Page *)left;
    hdr.level = level;
    records[0].key = key;
    records[0].ptr = (char *)right;
    records[1].ptr = NULL;

    hdr.last_index = 0;
    hdr.page_type = (uint8_t)page_type_;
    if (page_type_ == PageType::DRAM_INDEXTREE_PAGE || page_type_ == PageType::INDEXTREE_LAST_LEVEL_PAGE) {
#ifndef USE_SPINLOCK
      hdr.mtx = new std::mutex();
#else
      hdr.mtx = new Spinlock();
#endif
    }

    if (page_is_inpmem()) {
      clflush((char *)this, sizeof(Page));
    }
  }

  void *operator new(size_t size, bool inpmem);

  bool page_is_inpmem() {
    return hdr.page_type == PageType::NVM_SUBTREE_PAGE;
  }

  uint32_t GetLevel() {
    return hdr.level;
  }

  // 返回条目数量
  int count();

  bool remove_key(entry_key_t key);

  // remove a entry with no rebalance, need to acquire the lock
  bool remove(SubTree *bt, entry_key_t key, bool only_rebalance = false, bool with_lock = true);

  bool remove(IndexTree *bt, entry_key_t key, bool only_rebalance = false, bool with_lock = true);

  bool remove_rebalancing(SubTree *bt, entry_key_t key, bool only_rebalance = false, bool with_lock = true);

  void insert_key(entry_key_t key, char *ptr, int *num_entries, bool flush = true, bool update_last_index = true);

  // Insert a new key - FAST and FAIR
  Page *store(SubTree *bt, char *left, entry_key_t key, char *right, bool flush,
              bool with_lock, Page *invalid_sibling = NULL);

  Page * store(IndexTree *bt, char *left, entry_key_t key, char *right, bool flush,
              bool with_lock, Page *invalid_sibling = NULL);

  // Search keys with linear search
  void linear_search_range(entry_key_t min, entry_key_t max, unsigned long *buf);

  char *linear_search(entry_key_t key);

  // print a node
  void print();

  void printAll();
  
  void CalculateSapce(uint64_t &space);

  // Search keys with linear search
  void linear_search_range(entry_key_t min, entry_key_t max, 
      std::vector<std::pair<uint64_t, uint64_t>> &result, int &size);

  void linear_search_range(entry_key_t min, 
        entry_key_t max, void **values, int &size);
};

  void * Page::operator new(size_t size, bool inpmem) {
    void *ret = NULL;
    // ret = malloc(size);
    if (inpmem)
      alloc_memalign(&ret, 64, size);
    else {
      // int res = posix_memalign(&ret, 8, size);
      ret = malloc(size);
      assert(ret);
      // p_assert(res == 0, "alloc fail, res: %d, size: %ld", res, size);
    }
    return ret;
  }

  int Page::count() {
    uint8_t previous_switch_counter;
    int count = 0;
    do {
      previous_switch_counter = hdr.switch_counter;
      count = hdr.last_index + 1;

      while (count >= 0 && records[count].ptr != NULL) {
        if (IS_FORWARD(previous_switch_counter))
          ++count;
        else
          --count;
      }

      if (count < 0) {
        count = 0;
        while (records[count].ptr != NULL) {
          ++count;
        }
      }

    } while (previous_switch_counter != hdr.switch_counter); // 保证这个过程中没有新的update或insert操作

    return count;
  }

  bool Page::remove_key(entry_key_t key) {
    // Set the switch_counter
    if (IS_FORWARD(hdr.switch_counter))
      ++hdr.switch_counter;

    bool shift = false;
    int i;
    // 从左到右遍历搜索key
    for (i = 0; records[i].ptr != NULL; ++i) {
        // zzy add
      // NVM::const_stat.AddCompare();
      if (!shift && records[i].key == key) {
        records[i].ptr =
            (i == 0) ? (char *)hdr.leftmost_ptr : records[i - 1].ptr;
        shift = true;
      }
      // shift
      if (shift) {
        records[i].key = records[i + 1].key;
        records[i].ptr = records[i + 1].ptr;

        // flush
        if(page_is_inpmem()) {
          uint64_t records_ptr = (uint64_t)(&records[i]);
          int remainder = records_ptr % CACHE_LINE_SIZE;
          bool do_flush =
              (remainder == 0) ||
              ((((int)(remainder + sizeof(entry)) / CACHE_LINE_SIZE) == 1) &&
              ((remainder + sizeof(entry)) % CACHE_LINE_SIZE) != 0);
          if (do_flush) {
              clflush((char *)records_ptr, CACHE_LINE_SIZE);
          }
        }
      }
    }

    if (shift) {
      --hdr.last_index;
    //   zzy add
      if (page_is_inpmem()) {
        clflush((char *)&(hdr.last_index), sizeof(int8_t));
      }
    }
    return shift;
  }

  bool Page::remove(SubTree *bt, entry_key_t key, bool only_rebalance, bool with_lock) {
    if (with_lock) {
      hdr.mtx->lock();
    }

    bool ret = remove_key(key);

    if (with_lock) {
      hdr.mtx->unlock();
    }
    
    return ret;
  }

  // remove a entry with no rebalance, need to acquire the lock
  bool Page::remove(IndexTree *bt, entry_key_t key, bool only_rebalance, bool with_lock) {
    p_assert(false, "not implement yet");
    hdr.mtx->lock();

    bool ret = remove_key(key);

    hdr.mtx->unlock();

    return ret;
  }

  void Page::insert_key(entry_key_t key, char *ptr, int *num_entries, bool flush, bool update_last_index) {
    // update switch_counter
    // 在shift前,如果上次的操作是删除操作，则先更新switch_counter
    if (!IS_FORWARD(hdr.switch_counter))
      ++hdr.switch_counter;

    // FAST
    if (*num_entries == 0) { // this Page is empty
      entry *new_entry = (entry *)&records[0];
      entry *array_end = (entry *)&records[1];
      new_entry->key = (entry_key_t)key;
      new_entry->ptr = (char *)ptr;

      array_end->ptr = (char *)NULL;

      if (flush) {
        if (page_is_inpmem()) {
          clflush((char *)this, CACHE_LINE_SIZE);
        }
      }
    } else {
      int i = *num_entries - 1, inserted = 0, to_flush_cnt = 0;
      records[*num_entries + 1].ptr = records[*num_entries].ptr; 
      if (flush) {
        if ((uint64_t) & (records[*num_entries + 1].ptr) % CACHE_LINE_SIZE == 0) {
          if (page_is_inpmem()) {
            clflush((char *)&(records[*num_entries + 1].ptr), sizeof(char *));
          }
        }
      }

      // FAST
      for (i = *num_entries - 1; i >= 0; i--) {
        //   zzy add
          // NVM::const_stat.AddCompare();
        if (key < records[i].key) {
          records[i + 1].ptr = records[i].ptr;
          records[i + 1].key = records[i].key;

          if (flush) {
            uint64_t records_ptr = (uint64_t)(&records[i + 1]);

            int remainder = records_ptr % CACHE_LINE_SIZE;
            bool do_flush =
                (remainder == 0) ||
                ((((int)(remainder + sizeof(entry)) / CACHE_LINE_SIZE) == 1) &&
                 ((remainder + sizeof(entry)) % CACHE_LINE_SIZE) != 0);
            if (do_flush) {
              if (page_is_inpmem()) {
                clflush((char *)records_ptr, CACHE_LINE_SIZE);
              }
              to_flush_cnt = 0;
            } else
              ++to_flush_cnt;
          }
        } else {
          records[i + 1].ptr = records[i].ptr;
          records[i + 1].key = key;
          records[i + 1].ptr = ptr;

          if (flush) {
            if (page_is_inpmem()) {
              clflush((char *)&records[i + 1], sizeof(entry));
            }
          }
          inserted = 1;
          break;
        }
      }
      if (inserted == 0) {
        records[0].ptr = (char *)hdr.leftmost_ptr;
        records[0].key = key;
        records[0].ptr = ptr;
        if (flush) {
          if (page_is_inpmem()) {
            clflush((char *)&records[0], sizeof(entry));
          }
        }
      }
    }

    if (update_last_index) {
      hdr.last_index = (int8_t)(*num_entries);
      // zzy add
      if (page_is_inpmem()) {
        clflush((char *)&(hdr.last_index), sizeof(int8_t));
      }
    }
    ++(*num_entries);
  }

  
//   /*
//    * Although we implemented the rebalancing of B+-Tree, it is currently blocked
//    * for the performance. Please refer to the follow. Chi, P., Lee, W. C., &
//    * Xie, Y. (2014, August). Making B+-tree efficient in PCM-based main memory.
//    * In Proceedings of the 2014 international symposium on Low power electronics
//    * and design (pp. 69-74). ACM.
//    */
//   bool Page::remove_rebalancing(SubTree *bt, entry_key_t key, bool only_rebalance, bool with_lock) {
//     if (with_lock) {
//       hdr.mtx->lock();
//     }
//     if (hdr.is_deleted) {
//       if (with_lock) {
//         hdr.mtx->unlock();
//       }
//       return false;
//     }

//     if (!only_rebalance) {
//       int num_entries_before = count(); 

//       // This node is root
//       if (this == (Page *)bt->root) {
//         if (hdr.level > 0) { // 根节点非叶子节点
//           if (num_entries_before == 1 && !hdr.right_sibling_ptr) {
//             bt->root = (char *)hdr.leftmost_ptr;
//             clflush((char *)&(bt->root), sizeof(char *));

//             hdr.is_deleted = 1;
//           }
//         }

//         // Remove the key from this node
//         bool ret = remove_key(key);

//         if (with_lock) {
//           hdr.mtx->unlock();
//         }
//         return true;
//       }

//       bool should_rebalance = true;
//       // check the node utilization
//       if (num_entries_before - 1 >= (int)((cardinality - 1) * 0.5)) {
//         should_rebalance = false;
//       }

//       // Remove the key from this node
//       bool ret = remove_key(key);

//       if (!should_rebalance) {
//         if (with_lock) {
//           hdr.mtx->unlock();
//         }
//         return (hdr.leftmost_ptr == NULL) ? ret : true;
//       }
//     }

//     // Remove a key from the parent node
//     entry_key_t deleted_key_from_parent = 0;
//     bool is_leftmost_node = false;
//     Page *left_sibling;
//     bt->btree_delete_internal(key, (char *)this, hdr.level + 1,
//                               &deleted_key_from_parent, &is_leftmost_node,
//                               &left_sibling);

//     if (is_leftmost_node) {
//       if (with_lock) {
//         hdr.mtx->unlock();
//       }

//       if (!with_lock) {
//         hdr.right_sibling_ptr->hdr.mtx->lock();
//       }
//       hdr.right_sibling_ptr->remove(bt, hdr.right_sibling_ptr->records[0].key, true,
//                               with_lock);
//       if (!with_lock) {
//         hdr.right_sibling_ptr->hdr.mtx->unlock();
//       }
//       return true;
//     }

//     if (with_lock) {
//       left_sibling->hdr.mtx->lock();
//     }

//     while (left_sibling->hdr.right_sibling_ptr != this) {
//       if (with_lock) {
//         Page *t = left_sibling->hdr.right_sibling_ptr;
//         left_sibling->hdr.mtx->unlock();
//         left_sibling = t;
//         left_sibling->hdr.mtx->lock();
//       } else
//         left_sibling = left_sibling->hdr.right_sibling_ptr;
//     }

//     int num_entries = count();
//     int left_num_entries = left_sibling->count();

//     // Merge or Redistribution
//     int total_num_entries = num_entries + left_num_entries;
//     if (hdr.leftmost_ptr)
//       ++total_num_entries;

//     entry_key_t parent_key;

//     if (total_num_entries > cardinality - 1) { // Redistribution
//       int m = (int)ceil(total_num_entries / 2);

//       if (num_entries < left_num_entries) { // left -> right
//         if (hdr.leftmost_ptr == nullptr) {
//           for (int i = left_num_entries - 1; i >= m; i--) {
//             insert_key(left_sibling->records[i].key,
//                        left_sibling->records[i].ptr, &num_entries);
//           }

//           left_sibling->records[m].ptr = nullptr;
//           clflush((char *)&(left_sibling->records[m].ptr), sizeof(char *));

//           left_sibling->hdr.last_index = (int8_t)(m - 1);
//           clflush((char *)&(left_sibling->hdr.last_index), sizeof(int8_t));

//           parent_key = records[0].key;
//         } else {
//           insert_key(deleted_key_from_parent, (char *)hdr.leftmost_ptr,
//                      &num_entries);

//           for (int i = left_num_entries - 1; i > m; i--) {
//             insert_key(left_sibling->records[i].key,
//                        left_sibling->records[i].ptr, &num_entries);
//           }

//           parent_key = left_sibling->records[m].key;

//           hdr.leftmost_ptr = (Page *)left_sibling->records[m].ptr;
//           clflush((char *)&(hdr.leftmost_ptr), sizeof(Page *));

//           left_sibling->records[m].ptr = nullptr;
//           clflush((char *)&(left_sibling->records[m].ptr), sizeof(char *));

//           left_sibling->hdr.last_index = (int8_t)(m - 1);
//           clflush((char *)&(left_sibling->hdr.last_index), sizeof(int8_t));
//         }

//         if (left_sibling == ((Page *)bt->root)) {
//           Page *new_root =
//               new Page(left_sibling, parent_key, this, hdr.level + 1);
//           bt->setNewRoot((char *)new_root);
//         } else {
//           bt->btree_insert_internal((char *)left_sibling, parent_key,
//                                     (char *)this, hdr.level + 1);
//         }
//       } else { // from leftmost case
//         hdr.is_deleted = 1;
//         clflush((char *)&(hdr.is_deleted), sizeof(uint8_t));

//         Page *new_sibling = new Page(hdr.level);
//         new_sibling->hdr.mtx->lock();
//         new_sibling->hdr.right_sibling_ptr = hdr.right_sibling_ptr;

//         int num_dist_entries = num_entries - m;
//         int new_sibling_cnt = 0;

//         if (hdr.leftmost_ptr == nullptr) {
//           for (int i = 0; i < num_dist_entries; i++) {
//             left_sibling->insert_key(records[i].key, records[i].ptr,
//                                      &left_num_entries);
//           }

//           for (int i = num_dist_entries; records[i].ptr != NULL; i++) {
//             new_sibling->insert_key(records[i].key, records[i].ptr,
//                                     &new_sibling_cnt, false);
//           }

//           clflush((char *)(new_sibling), sizeof(Page));

//           left_sibling->hdr.right_sibling_ptr = new_sibling;
//           clflush((char *)&(left_sibling->hdr.right_sibling_ptr), sizeof(Page *));

//           parent_key = new_sibling->records[0].key;
//         } else {
//           left_sibling->insert_key(deleted_key_from_parent,
//                                    (char *)hdr.leftmost_ptr, &left_num_entries);

//           for (int i = 0; i < num_dist_entries - 1; i++) {
//             left_sibling->insert_key(records[i].key, records[i].ptr,
//                                      &left_num_entries);
//           }

//           parent_key = records[num_dist_entries - 1].key;

//           new_sibling->hdr.leftmost_ptr =
//               (Page *)records[num_dist_entries - 1].ptr;
//           for (int i = num_dist_entries; records[i].ptr != NULL; i++) {
//             new_sibling->insert_key(records[i].key, records[i].ptr,
//                                     &new_sibling_cnt, false);
//           }
//           clflush((char *)(new_sibling), sizeof(Page));

//           left_sibling->hdr.right_sibling_ptr = new_sibling;
//           clflush((char *)&(left_sibling->hdr.right_sibling_ptr), sizeof(Page *));
//         }

//         if (left_sibling == ((Page *)bt->root)) {
//           Page *new_root =
//               new Page(left_sibling, parent_key, new_sibling, hdr.level + 1);
//           bt->setNewRoot((char *)new_root);
//         } else {
//           bt->btree_insert_internal((char *)left_sibling, parent_key,
//                                     (char *)new_sibling, hdr.level + 1);
//         }

//         new_sibling->hdr.mtx->unlock();
//       }
//     } else {
//       hdr.is_deleted = 1;
//       clflush((char *)&(hdr.is_deleted), sizeof(uint8_t));

//       if (hdr.leftmost_ptr)
//         left_sibling->insert_key(deleted_key_from_parent,
//                                  (char *)hdr.leftmost_ptr, &left_num_entries);

//       for (int i = 0; records[i].ptr != NULL; ++i) {
//         left_sibling->insert_key(records[i].key, records[i].ptr,
//                                  &left_num_entries);
//       }

//       left_sibling->hdr.right_sibling_ptr = hdr.right_sibling_ptr;
//       clflush((char *)&(left_sibling->hdr.right_sibling_ptr), sizeof(Page *));
//     }

//     if (with_lock) {
//       left_sibling->hdr.mtx->unlock();
//       hdr.mtx->unlock();
//     }

//     return true;
//   }

  Page * Page::store(SubTree *bt, char *left, entry_key_t key, char *right, bool flush,
              bool with_lock, Page *invalid_sibling) {
    if (with_lock) {
      hdr.mtx->lock(); // Lock the write lock
    }
    if (hdr.is_deleted) {
      if (with_lock) {
        hdr.mtx->unlock();
      }

      return NULL;
    }

    // If this node has a sibling node,
    if (hdr.right_sibling_ptr && (hdr.right_sibling_ptr != invalid_sibling)) {
      // Compare this key with the first key of the sibling
      if (key > hdr.right_sibling_ptr->records[0].key) { 
        // 如果兄弟节点的最小key大于要插入的key,
        // 则在兄弟节点执行插入，对应论文Fair算法中更新父节点时出现并发Insert的情况
        // zzy add
        // NVM::const_stat.AddCompare();
        if (with_lock) {
          hdr.mtx->unlock(); // Unlock the write lock
        }
        // return hdr.right_sibling_ptr->store(bt, NULL, key, right, true, with_lock,
        //                               invalid_sibling);
        return NULL;
      }
    }

    int num_entries = count();

    // FAST
    if (num_entries < cardinality - 1) { // no need to split
      insert_key(key, right, &num_entries, flush);

      if (with_lock) {
        hdr.mtx->unlock(); // Unlock the write lock
      }

      return this;
    } else { // FAIR, need to split
      // overflow
      // create a new node
      Page *sibling = nullptr;
      if (page_is_inpmem()) {
        sibling = new(true) Page(PageType::NVM_SUBTREE_PAGE, hdr.level);
      } else {
        sibling = new(false) Page((PageType)this->hdr.page_type, hdr.level);
      }
      
      int m = (int)ceil(num_entries / 2);
      entry_key_t split_key = records[m].key;

      // migrate half of keys into the sibling
      int sibling_cnt = 0;
      if (hdr.leftmost_ptr == NULL) { // leaf node
        for (int i = m; i < num_entries; ++i) {
          sibling->insert_key(records[i].key, records[i].ptr, &sibling_cnt,
                              false);
        }
      } else { // internal node
        for (int i = m + 1; i < num_entries; ++i) { // XXX: 这里i从m+1似乎并不平均，假设num_entries=3，则实际没移动
          sibling->insert_key(records[i].key, records[i].ptr, &sibling_cnt,
                              false);
        }
        sibling->hdr.leftmost_ptr = (Page *)records[m].ptr;
      }

      sibling->hdr.right_sibling_ptr = hdr.right_sibling_ptr;
      if (page_is_inpmem()) {
        clflush((char *)sibling, sizeof(Page));
      }
      
      hdr.right_sibling_ptr = sibling;
      if (page_is_inpmem()) {
        clflush((char *)&hdr, sizeof(hdr));
      }

      // set to NULL
      if (IS_FORWARD(hdr.switch_counter)) // XXX: 注意insert操作在设置完sibling指针即增加switch_counter
        hdr.switch_counter += 2;
      else
        ++hdr.switch_counter;
      records[m].ptr = NULL;

      if (page_is_inpmem()) {
        clflush((char *)&records[m], sizeof(entry));
      }

      hdr.last_index = (int8_t)(m - 1);
      if (page_is_inpmem()) {
        clflush((char *)&(hdr.last_index), sizeof(int8_t));
      }
      
      num_entries = hdr.last_index + 1;

      Page *ret;

      // insert the key
      // zzy add
      // NVM::const_stat.AddCompare();
      if (key < split_key) {
        insert_key(key, right, &num_entries);
        ret = this;
      } else {
        sibling->insert_key(key, right, &sibling_cnt);
        ret = sibling;
      }

      if (hdr.page_type == PageType::NVM_SUBTREE_PAGE) {
        // Set a new root or insert the split key to the parent
        if (bt->sub_root_ == (char *)this) { // only one node can update the root ptr
          // zzy add
          // 如果子树高度超过MAX_SUBTREE_HEIGHT
            // 如果内存还没创建index_tree则创建
            // 否则分裂subtree，插入index_tree， index_tree最后一层存的是SubTree的指针
          if (hdr.level + 1 >= MAX_SUBTREE_HEIGHT) {
            SubTree * sibling_subtree = new SubTree(sibling, SubTreeStatus::IN_NVM);
            sibling_subtree->left_sibling_subtree_ = bt;
            sibling_subtree->right_sibling_subtree_ = bt->right_sibling_subtree_;
            bt->right_sibling_subtree_ = sibling_subtree;
            if (sibling_subtree->right_sibling_subtree_)
              sibling_subtree->right_sibling_subtree_->left_sibling_subtree_ = sibling_subtree;
            for (int i = 0; i < numa_node_num; i++) {
              bt->hotness_[i] /= 2;
              sibling_subtree->hotness_[i] = bt->hotness_[i];
            }
            
            for (int nn = 0; nn < numa_node_num; nn++) {
              bt->read_times_[nn] /= 2;
              sibling_subtree->read_times_[nn] = bt->read_times_[nn];

              bt->write_times_[nn] /= 2;
              sibling_subtree->write_times_[nn] = bt->write_times_[nn];
            }

            if (unlikely(!index_tree_root->has_hasindextree())) {
              index_tree_root->indextree_lock();
              if (index_tree_root->has_hasindextree()) {
                index_tree_root->indextree_unlock();
              }
            }

            if (!(index_tree_root->has_hasindextree())) {
              Page *index_root_page = new(false) Page(PageType::INDEXTREE_LAST_LEVEL_PAGE, bt, split_key, sibling_subtree, hdr.level + 1);
              IndexTree *indextree_ = new IndexTree(index_root_page, hdr.level + 1);
              index_tree_root->set_indextree(indextree_);
              index_tree_root->indextree_unlock();
              if (with_lock) {
                hdr.mtx->unlock(); // Unlock the write lock
              }
              std::cout << "generate index tree: " << index_tree_root << std::endl << std::flush;
            } else {
              if (with_lock) {
                hdr.mtx->unlock(); // Unlock the write lock
              }
              // 调整父节点
              IndexTree *indextree_ = (IndexTree *)(index_tree_root->tree_);
              indextree_->btree_insert_internal(NULL, split_key, (char *)sibling_subtree, hdr.level + 1);
            }
            statis_->insert_subtree(sibling_subtree);
          } else {
            Page *new_root = new(true) Page(PageType::NVM_SUBTREE_PAGE, (Page *)this, split_key, sibling, hdr.level + 1);
            bt->setNewRoot((char *)new_root);

            if (with_lock) {
              hdr.mtx->unlock(); // Unlock the write lock
            }
          }
          // Page *new_root =
          //     new Page((Page *)this, split_key, sibling, hdr.level + 1);
          // bt->setNewRoot((char *)new_root);

          // if (with_lock) {
          //   hdr.mtx->unlock(); // Unlock the write lock
          // }
        } else {
          if (with_lock) {
            hdr.mtx->unlock(); // Unlock the write lock
          }
          // 调整父节点
          bt->btree_insert_internal(NULL, split_key, (char *)sibling, hdr.level + 1);
        }
        
      } else if (hdr.page_type == PageType::DRAM_CACHETREE_PAGE) {
        // Set a new root or insert the split key to the parent
        if (bt->sub_root_ == (char *)this) { // only one node can update the root ptr
          // zzy add
          // 如果子树高度超过MAX_SUBTREE_HEIGHT
            // 如果内存还没创建index_tree则创建
            // 否则分裂subtree，插入index_tree， index_tree最后一层存的是SubTree的指针
          if (hdr.level + 1 >= MAX_SUBTREE_HEIGHT) {
            SubTree * sibling_subtree = new SubTree(sibling, SubTreeStatus::IN_DRAM);
            sibling_subtree->treelog_ = new TreeLog();
            sibling_subtree->left_sibling_subtree_ = bt;
            sibling_subtree->right_sibling_subtree_ = bt->right_sibling_subtree_;
            bt->right_sibling_subtree_ = sibling_subtree;
            if (sibling_subtree->right_sibling_subtree_)
              sibling_subtree->right_sibling_subtree_->left_sibling_subtree_ = sibling_subtree;
            for (int i = 0; i < numa_node_num; i++) {
              bt->hotness_[i] /= 2;
              sibling_subtree->hotness_[i] = bt->hotness_[i];
            }
            
            for (int nn = 0; nn < numa_node_num; nn++) {
              bt->read_times_[nn] /= 2;
              sibling_subtree->read_times_[nn] = bt->read_times_[nn];

              bt->write_times_[nn] /= 2;
              sibling_subtree->write_times_[nn] = bt->write_times_[nn];
            }

            if (unlikely(!index_tree_root->has_hasindextree())) {
              index_tree_root->indextree_lock();
              if (index_tree_root->has_hasindextree()) {
                index_tree_root->indextree_unlock();
              }
            }

            if (!(index_tree_root->has_hasindextree())) {
              Page *index_root_page = new(false) Page(PageType::INDEXTREE_LAST_LEVEL_PAGE, bt, split_key, sibling_subtree, hdr.level + 1);
              IndexTree *indextree_ = new IndexTree(index_root_page, hdr.level + 1);
              index_tree_root->set_indextree(indextree_);
              index_tree_root->indextree_unlock();
              if (with_lock) {
                hdr.mtx->unlock(); // Unlock the write lock
              }
              std::cout << "generate index tree: " << index_tree_root << std::endl << std::flush;
            } else {
              if (with_lock) {
                hdr.mtx->unlock(); // Unlock the write lock
              }
              // 调整父节点
              IndexTree *indextree_ = (IndexTree *)(index_tree_root->tree_);
              indextree_->btree_insert_internal(NULL, split_key, (char *)sibling_subtree, hdr.level + 1);
            }
            statis_->insert_subtree(sibling_subtree);
          } else {
            Page *new_root = new(false) Page(PageType::DRAM_CACHETREE_PAGE, (Page *)this, split_key, sibling, hdr.level + 1);
            bt->setNewRoot((char *)new_root);

            if (with_lock) {
              hdr.mtx->unlock(); // Unlock the write lock
            }
          }
          // Page *new_root =
          //     new Page((Page *)this, split_key, sibling, hdr.level + 1);
          // bt->setNewRoot((char *)new_root);

          // if (with_lock) {
          //   hdr.mtx->unlock(); // Unlock the write lock
          // }
        } else {
          if (with_lock) {
            hdr.mtx->unlock(); // Unlock the write lock
          }
          // 调整父节点
          bt->btree_insert_internal(NULL, split_key, (char *)sibling, hdr.level + 1);
        }
      } else {
        p_assert(false, "should never be here");
      }
      return ret;
    }
  }

  // Insert a new key - FAST and FAIR
  Page * Page::store(IndexTree *bt, char *left, entry_key_t key, char *right, bool flush,
              bool with_lock, Page *invalid_sibling) {
    if (with_lock) {
      hdr.mtx->lock(); // Lock the write lock
    }
    if (hdr.is_deleted) {
      if (with_lock) {
        hdr.mtx->unlock();
      }

      return NULL;
    }

    // If this node has a sibling node,
    if (hdr.right_sibling_ptr && (hdr.right_sibling_ptr != invalid_sibling)) {
      // Compare this key with the first key of the sibling
      if (key > hdr.right_sibling_ptr->records[0].key) { 
        // 如果兄弟节点的最小key大于要插入的key,
        // 则在兄弟节点执行插入，对应论文Fair算法中更新父节点时出现并发Insert的情况
        // zzy add
        // NVM::const_stat.AddCompare();
        if (with_lock) {
          hdr.mtx->unlock(); // Unlock the write lock
        }
        return hdr.right_sibling_ptr->store(bt, NULL, key, right, true, with_lock,
                                      invalid_sibling);
      }
    }

    int num_entries = count();

    // FAST
    if (num_entries < cardinality - 1) { // no need to split
      insert_key(key, right, &num_entries, flush);

      if (with_lock) {
        hdr.mtx->unlock(); // Unlock the write lock
      }

      return this;
    } else { // FAIR, need to split
      // overflow
      // create a new node
      Page *sibling = nullptr;
      p_assert(!page_is_inpmem(), "index tree must be in dram");
      sibling = new(false) Page((PageType)this->hdr.page_type, hdr.level);
      
      int m = (int)ceil(num_entries / 2);
      entry_key_t split_key = records[m].key;

      // migrate half of keys into the sibling
      int sibling_cnt = 0;
      if (hdr.leftmost_ptr == NULL) { // leaf node
        p_assert(false, "should not be here!");
        for (int i = m; i < num_entries; ++i) {
          sibling->insert_key(records[i].key, records[i].ptr, &sibling_cnt,
                              false);
        }
      } else { // internal node
        for (int i = m + 1; i < num_entries; ++i) { // XXX: 这里i从m+1似乎并不平均，假设num_entries=3，则实际没移动
          sibling->insert_key(records[i].key, records[i].ptr, &sibling_cnt,
                              false);
        }
        sibling->hdr.leftmost_ptr = (Page *)records[m].ptr;
      }

      sibling->hdr.right_sibling_ptr = hdr.right_sibling_ptr;
      
      hdr.right_sibling_ptr = sibling;

      // set to NULL
      if (IS_FORWARD(hdr.switch_counter)) // XXX: 注意insert操作在设置完sibling指针即增加switch_counter
        hdr.switch_counter += 2;
      else
        ++hdr.switch_counter;
      records[m].ptr = NULL;

      hdr.last_index = (int8_t)(m - 1);
      
      num_entries = hdr.last_index + 1;

      Page *ret;

      // insert the key
      // zzy add
      // NVM::const_stat.AddCompare();
      if (key < split_key) {
        insert_key(key, right, &num_entries);
        ret = this;
      } else {
        sibling->insert_key(key, right, &sibling_cnt);
        ret = sibling;
      }

      p_assert(hdr.page_type == PageType::DRAM_INDEXTREE_PAGE ||
        hdr.page_type == PageType::INDEXTREE_LAST_LEVEL_PAGE, "must be indextree_page");

      // Set a new root or insert the split key to the parent
      IndexTree *indextree_ = (IndexTree *)index_tree_root->tree_;
      if (indextree_->getRoot() == this) { // only one node can update the root ptr
        Page *index_root_page = new(false) Page(PageType::DRAM_INDEXTREE_PAGE, (Page *)this, split_key, sibling, hdr.level + 1);
        indextree_->setNewRoot((char *)index_root_page);
        // index_tree_root.
        if (with_lock) {
          hdr.mtx->unlock(); // Unlock the write lock
        }
      } else {
        if (with_lock) {
          hdr.mtx->unlock(); // Unlock the write lock
        }
        // 调整父节点
        p_assert(bt == indextree_, "should be equal");
        bt->btree_insert_internal(NULL, split_key, (char *)sibling, hdr.level + 1);
      }

      return ret;
    }
  }

  void Page::linear_search_range(entry_key_t min, entry_key_t max, unsigned long *buf) {
    int i, off = 0;
    uint8_t previous_switch_counter;
    Page *current = this;

    while (current) {
      int old_off = off;
      do {
        previous_switch_counter = current->hdr.switch_counter;
        off = old_off;

        entry_key_t tmp_key;
        char *tmp_ptr;

        if (IS_FORWARD(previous_switch_counter)) { // 
          if ((tmp_key = current->records[0].key) > min) {
            //   zzy add
            // NVM::const_stat.AddCompare();
            if (tmp_key < max) {
              if ((tmp_ptr = current->records[0].ptr) != NULL) {
                if (tmp_key == current->records[0].key) {
                  if (tmp_ptr) {
                    buf[off++] = (unsigned long)tmp_ptr;
                  }
                }
              }
            } else
              return;
          }

          for (i = 1; current->records[i].ptr != NULL; ++i) {
            if ((tmp_key = current->records[i].key) > min) {
                // zzy add
                // NVM::const_stat.AddCompare();
              if (tmp_key < max) {
                if ((tmp_ptr = current->records[i].ptr) !=
                    current->records[i - 1].ptr) {
                  if (tmp_key == current->records[i].key) {
                    if (tmp_ptr)
                      buf[off++] = (unsigned long)tmp_ptr;
                  }
                }
              } else
                return;
            }
          }
        } else {
          for (i = count() - 1; i > 0; --i) {
            if ((tmp_key = current->records[i].key) > min) {
                // zzy add
                // NVM::const_stat.AddCompare();
              if (tmp_key < max) {
                if ((tmp_ptr = current->records[i].ptr) !=
                    current->records[i - 1].ptr) {
                  if (tmp_key == current->records[i].key) {
                    if (tmp_ptr)
                      buf[off++] = (unsigned long)tmp_ptr;
                  }
                }
              } else
                return;
            }
          }

          if ((tmp_key = current->records[0].key) > min) {
            //   zzy add
            // NVM::const_stat.AddCompare();
            if (tmp_key < max) {
              if ((tmp_ptr = current->records[0].ptr) != NULL) {
                if (tmp_key == current->records[0].key) {
                  if (tmp_ptr) {
                    buf[off++] = (unsigned long)tmp_ptr;
                  }
                }
              }
            } else
              return;
          }
        }
      } while (previous_switch_counter != current->hdr.switch_counter); // 发生了修改操作，重试

      current = current->hdr.right_sibling_ptr; // 移动到下一个page
    }
  }

  char * Page::linear_search(entry_key_t key) {
    int i = 1;
    uint8_t previous_switch_counter;
    char *ret = NULL;
    char *t;
    entry_key_t k;

    if (hdr.leftmost_ptr == NULL) { // Search a leaf node
      do {
        previous_switch_counter = hdr.switch_counter;
        ret = NULL;

        // search from left ro right
        if (IS_FORWARD(previous_switch_counter)) {
            // zzy add
            // NVM::const_stat.AddCompare();
          if ((k = records[0].key) == key) {
            if ((t = records[0].ptr) != NULL) {
              if (k == records[0].key) {
                ret = t;
                continue;
              }
            }
          }

          for (i = 1; records[i].ptr != NULL; ++i) {
            //   zzy add
            // NVM::const_stat.AddCompare();
            if ((k = records[i].key) == key) {
              if (records[i - 1].ptr != (t = records[i].ptr)) {
                if (k == records[i].key) {
                  ret = t;
                  break;
                }
              }
            }
          }
        } else { // search from right to left
          for (i = count() - 1; i > 0; --i) {
            //   zzy add
            // NVM::const_stat.AddCompare();
            if ((k = records[i].key) == key) {
              if (records[i - 1].ptr != (t = records[i].ptr) && t) {
                if (k == records[i].key) {
                  ret = t;
                  break;
                }
              }
            }
          }

          if (!ret) {
            //   zzy add
            // NVM::const_stat.AddCompare();
            if ((k = records[0].key) == key) {
              if (NULL != (t = records[0].ptr) && t) {
                if (k == records[0].key) {
                  ret = t;
                  continue;
                }
              }
            }
          }
        }
      } while (hdr.switch_counter != previous_switch_counter);

      if (ret) {
        return ret;
      }

    //   zzy add
    // NVM::const_stat.AddCompare();

      if ((t = (char *)hdr.right_sibling_ptr) && key >= ((Page *)t)->records[0].key) {
        return t;
      }
        
      return NULL;
    } else { // internal node
      do {
        previous_switch_counter = hdr.switch_counter;
        ret = NULL;

        if (IS_FORWARD(previous_switch_counter)) {
            // zzy add
            // NVM::const_stat.AddCompare();
          if (key < (k = records[0].key)) {
            if ((t = (char *)hdr.leftmost_ptr) != records[0].ptr) {
              ret = t;
              continue;
            }
          }

          for (i = 1; records[i].ptr != NULL; ++i) {
            //   zzy add
            // NVM::const_stat.AddCompare();
            if (key < (k = records[i].key)) {
              if ((t = records[i - 1].ptr) != records[i].ptr) {
                ret = t;
                break;
              }
            }
          }

          if (!ret) {
            ret = records[i - 1].ptr;
            continue;
          }
        } else { // search from right to left
          for (i = count() - 1; i >= 0; --i) {
            //   zzy add
            // NVM::const_stat.AddCompare();
            if (key >= (k = records[i].key)) {
              if (i == 0) {
                if ((char *)hdr.leftmost_ptr != (t = records[i].ptr)) {
                  ret = t;
                  break;
                }
              } else {
                if (records[i - 1].ptr != (t = records[i].ptr)) {
                  ret = t;
                  break;
                }
              }
            }
          }
        }
      } while (hdr.switch_counter != previous_switch_counter);

      if ((t = (char *)hdr.right_sibling_ptr) != NULL) {
        //   zzy add
        // NVM::const_stat.AddCompare();
        if (key >= ((Page *)t)->records[0].key) {
          return t;
        }
      }

      if (ret) {
        return ret;
      } else
        return (char *)hdr.leftmost_ptr;
    }

    return NULL;
  }

  void Page::print() {
    if (hdr.leftmost_ptr == NULL)
      printf("[%d] leaf %p \n", this->hdr.level, this);
    else
      printf("[%d] internal %p \n", this->hdr.level, this);
    printf("last_index: %d\n", hdr.last_index);
    printf("switch_counter: %d\n", hdr.switch_counter);
    printf("search direction: ");
    if (IS_FORWARD(hdr.switch_counter))
      printf("->\n");
    else
      printf("<-\n");

    if (hdr.leftmost_ptr != NULL)
      printf("%p ", hdr.leftmost_ptr);

    for (int i = 0; records[i].ptr != NULL; ++i)
      printf("%ld,%p ", records[i].key, records[i].ptr);

    printf("%p ", hdr.right_sibling_ptr);

    printf("\n");
  }

  void Page::printAll() {
    if (hdr.leftmost_ptr == NULL) {
      printf("printing leaf node: ");
      print();
    } else {
      printf("printing internal node: ");
      print();
      ((Page *)hdr.leftmost_ptr)->printAll();
      for (int i = 0; records[i].ptr != NULL; ++i) {
        ((Page *)records[i].ptr)->printAll();
      }
    }
  }

  void Page::CalculateSapce(uint64_t &space) {
      if(hdr.leftmost_ptr==NULL) {
          space += PAGESIZE;
      } else {
          space += PAGESIZE;
          ((Page*) hdr.leftmost_ptr)->CalculateSapce(space);
          for(int i=0;records[i].ptr != NULL;++i){
              ((Page*) records[i].ptr)->CalculateSapce(space);
          }
      }
  }

  void Page::linear_search_range(entry_key_t min, entry_key_t max, 
      std::vector<std::pair<uint64_t, uint64_t>> &result, int &size) {
    int i, off = 0;
    uint8_t previous_switch_counter;
    Page *current = this;

    while (current) {
      int old_off = off;
      do {
        previous_switch_counter = current->hdr.switch_counter;
        off = old_off;

        entry_key_t tmp_key;
        char *tmp_ptr;

        if (IS_FORWARD(previous_switch_counter)) { // 
          if ((tmp_key = current->records[0].key) > min) {
            //   zzy add
            // NVM::const_stat.AddCompare();
            if (tmp_key < max) {
              if ((tmp_ptr = current->records[0].ptr) != NULL) {
                if (tmp_key == current->records[0].key) {
                  if (tmp_ptr) {
                    // zzy add
                    // buf[off++] = (unsigned long)tmp_ptr;
                    result.push_back({tmp_key, (uint64_t)tmp_ptr});
                    off++;
                    if(off >= size) {
                      return;
                    }
                  }
                }
              }
            } else {
              // zzy add
              size = off;
              return;
            }
          }

          for (i = 1; current->records[i].ptr != NULL; ++i) {
            if ((tmp_key = current->records[i].key) > min) {
              // zzy add
              // NVM::const_stat.AddCompare();
              if (tmp_key < max) {
                if ((tmp_ptr = current->records[i].ptr) !=
                    current->records[i - 1].ptr) {
                  if (tmp_key == current->records[i].key) {
                    if (tmp_ptr) {
                      // buf[off++] = (unsigned long)tmp_ptr;
                      result.push_back({tmp_key, (uint64_t)tmp_ptr});
                      off++;
                      if(off >= size) {
                        return;
                      }
                    }
                      
                  }
                }
              } else {
                // zzy add
                size = off;
                return;
              } 
            }
          }
        } else {
          for (i = count() - 1; i > 0; --i) {
            if ((tmp_key = current->records[i].key) > min) {
              // zzy add
              // NVM::const_stat.AddCompare();
              if (tmp_key < max) {
                if ((tmp_ptr = current->records[i].ptr) !=
                    current->records[i - 1].ptr) {
                  if (tmp_key == current->records[i].key) {
                    if (tmp_ptr) {
                      // buf[off++] = (unsigned long)tmp_ptr;
                      // zzy add
                      result.push_back({tmp_key, (uint64_t)tmp_ptr});
                      off++;
                      if(off >= size) {
                        return;
                      }
                    }
                  }
                }
              } else {
                size = off;
                return;
              }
            }
          }

          if ((tmp_key = current->records[0].key) > min) {
            //   zzy add
            // NVM::const_stat.AddCompare();
            if (tmp_key < max) {
              if ((tmp_ptr = current->records[0].ptr) != NULL) {
                if (tmp_key == current->records[0].key) {
                  if (tmp_ptr) {
                    // buf[off++] = (unsigned long)tmp_ptr;
                    // zzy add
                    result.push_back({tmp_key, (uint64_t)tmp_ptr});
                    off++;
                    if(off >= size) {
                      return;
                    }
                  }
                }
              }
            } else {
              size = off;
              return;
            }
          }
        }
      } while (previous_switch_counter != current->hdr.switch_counter); // 发生了修改操作，重试

      current = current->hdr.right_sibling_ptr; // 移动到下一个page
    }
    size = off;
  }

  void Page::linear_search_range(entry_key_t min, 
        entry_key_t max, void **values, int &size) {
    int i, off = 0;
    uint8_t previous_switch_counter;
    Page *current = this;

    while (current) {
      int old_off = off;
      do {
        previous_switch_counter = current->hdr.switch_counter;
        off = old_off;

        entry_key_t tmp_key;
        char *tmp_ptr;

        if (IS_FORWARD(previous_switch_counter)) { // 
          if ((tmp_key = current->records[0].key) > min) {
            //   zzy add
            // NVM::const_stat.AddCompare();
            if (tmp_key < max) {
              if ((tmp_ptr = current->records[0].ptr) != NULL) {
                if (tmp_key == current->records[0].key) {
                  if (tmp_ptr) {
                    // zzy add
                    // buf[off++] = (unsigned long)tmp_ptr;
                    values[off++] = tmp_ptr;
                    if(off >= size) {
                      return;
                    }
                  }
                }
              }
            } else {
              // zzy add
              size = off;
              return;
            }
          }

          for (i = 1; current->records[i].ptr != NULL; ++i) {
            if ((tmp_key = current->records[i].key) > min) {
              // zzy add
              // NVM::const_stat.AddCompare();
              if (tmp_key < max) {
                if ((tmp_ptr = current->records[i].ptr) !=
                    current->records[i - 1].ptr) {
                  if (tmp_key == current->records[i].key) {
                    if (tmp_ptr) {
                      // buf[off++] = (unsigned long)tmp_ptr;
                      values[off++] = tmp_ptr;
                      if(off >= size) {
                        return;
                      }
                    }
                      
                  }
                }
              } else {
                // zzy add
                size = off;
                return;
              } 
            }
          }
        } else {
          for (i = count() - 1; i > 0; --i) {
            if ((tmp_key = current->records[i].key) > min) {
              // zzy add
              // NVM::const_stat.AddCompare();
              if (tmp_key < max) {
                if ((tmp_ptr = current->records[i].ptr) !=
                    current->records[i - 1].ptr) {
                  if (tmp_key == current->records[i].key) {
                    if (tmp_ptr) {
                      // buf[off++] = (unsigned long)tmp_ptr;
                      // zzy add
                      values[off++] = tmp_ptr;
                      if(off >= size) {
                        return;
                      }
                    }
                  }
                }
              } else {
                size = off;
                return;
              }
            }
          }

          if ((tmp_key = current->records[0].key) > min) {
            //   zzy add
            // NVM::const_stat.AddCompare();
            if (tmp_key < max) {
              if ((tmp_ptr = current->records[0].ptr) != NULL) {
                if (tmp_key == current->records[0].key) {
                  if (tmp_ptr) {
                    // buf[off++] = (unsigned long)tmp_ptr;
                    // zzy add
                      values[off++] = tmp_ptr;
                      if(off >= size) {
                        return;
                      }
                  }
                }
              }
            } else {
              size = off;
              return;
            }
          }
        }
      } while (previous_switch_counter != current->hdr.switch_counter); // 发生了修改操作，重试

      current = current->hdr.right_sibling_ptr; // 移动到下一个page
    }
    size = off;
  }
}