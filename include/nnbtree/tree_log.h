#pragma once

#include "bitmap.h"
#include "nvm_alloc.h"
#include "numa_config.h"

#include <assert.h>
#include <atomic>

namespace nnbtree {

static inline void clflush(char *data, int len) {
    NVM::Mem_persist(data, len);
}

const uint64_t LogSize = (1UL << 18); // 每个子树日志256KB
const uint64_t LogNum = 65536; // 每个numa节点分配65536个Log

/**
 * @brief 每个日志的结构
 * 
 * +------------Log--------------+
 * +-num_logs(4B)-|----logs------+
 */

// 全局的log池，管理所有NUMA节点的log分配和释放
class TreeLogPool {
    public:
        TreeLogPool() {
            // 在每个NUMA节点上都分配
            memset(bitmap_, 0, sizeof(bitmap_));
            mapped_len_ = LogSize * LogNum;
            for (int i = 0; i < numa_node_num; i++) {
                bitmap_[i] = create_bitmap(LogNum);
            }
            assert(bitmap_);
            log_start_addr_ = index_pemm_alloc_log(mapped_len_);
        }

        ~TreeLogPool() {
            for (int i = 0; i < numa_node_num; i++) {
                pmem_unmap(log_start_addr_[i], mapped_len_);
            }

            free((void*)log_start_addr_);

            free(bitmap_);
        }

        char * foregroundthread_alloc_logfile() {
            bitmap * thread_bitmap = bitmap_[numa_map[my_thread_id]];

            if (bitmap_full(thread_bitmap)) {
                std::cout << "log pool is full" << std::endl;
                exit(-1);
            }

            int pos = get_free(thread_bitmap);

            char *log_addr = (char*)(log_start_addr_[numa_map[my_thread_id]]) + pos * LogSize;

            *((uint32_t *)log_addr) = 0;
            p_assert(*(uint32_t *)log_addr == 0, "init log metadata fail"); 

            return log_addr;
        }

        void foregroundthread_delete_logfile(char *log_addr) {
            bitmap * thread_bitmap = bitmap_[numa_map[my_thread_id]];
            int pos = (log_addr - (char *)log_start_addr_[numa_map[my_thread_id]]) / LogSize;
            put_back(thread_bitmap, pos);
        }

    private:
        void **log_start_addr_;
        bitmap *bitmap_[numa_node_num];
        size_t mapped_len_;
};

extern TreeLogPool *treelog_pool;

void init_treelogpool();

enum TreeLogType : uint8_t {
    INSERT,
    UPDATE,
    DELETE
};

struct TreeLogEntry {
    TreeLogType log_type_;
    uint64_t key_;
    uint64_t value_;
};

static_assert(sizeof(TreeLogEntry) == 24);

class TreeLog {
    public:
        TreeLog() :log_num_(0) {
            start_addr_ = treelog_pool->foregroundthread_alloc_logfile();
            cur_addr = start_addr_ + sizeof(uint32_t);
        }

        ~TreeLog() {
            treelog_pool->foregroundthread_delete_logfile(start_addr_);
        }

        void set_log_num(uint64_t log_num) {
            *((uint32_t *)start_addr_) = log_num;
            clflush(start_addr_, sizeof(uint32_t));
            log_num_ = log_num;
        }

        void write_log(TreeLogType log_type, uint64_t key, uint64_t val) {
            ((TreeLogEntry *)cur_addr)->log_type_ = log_type;
            ((TreeLogEntry *)cur_addr)->key_ = key;
            ((TreeLogEntry *)cur_addr)->value_ = val;
            clflush(cur_addr, sizeof(TreeLogEntry));
            log_num_++;
            *((uint32_t *)cur_addr) = log_num_;
            clflush(cur_addr, sizeof(TreeLogEntry));
        }

    private:
        char *start_addr_; // 从内存池分配的日志的起始位置，开头是元数据
        char *cur_addr; // 下一次要写log的addr
        uint32_t log_num_; // 日志数目
};


} // end namespace nnbtree