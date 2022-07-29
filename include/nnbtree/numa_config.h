#pragma once 

#include <libpmemobj++/detail/common.hpp>
#include <libpmemobj++/detail/template_helpers.hpp>
#include <libpmemobj++/experimental/v.hpp>
#include <libpmemobj++/make_persistent.hpp>
#include <libpmemobj++/make_persistent_atomic.hpp>
#include <libpmemobj++/mutex.hpp>
#include <libpmemobj++/p.hpp>
#include <libpmemobj++/persistent_ptr.hpp>
#include <libpmemobj++/transaction.hpp>

#include <pthread.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <sys/ioctl.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <unistd.h>

namespace nnbtree {

const int numa_max_node = 8; 

const int numa_node_num = 2; // current machine numa nodes

const int max_thread_num = 64; // 每个numa节点32threads

extern pmem::obj::pool_base pop_numa[nnbtree::numa_max_node];
extern int8_t numa_map[nnbtree::max_thread_num];
extern thread_local int my_thread_id;

void bindCore(uint16_t core);
void init_numa_pool();

void *index_pmem_alloc(size_t size);
void index_pmem_free(void *ptr);

// 在每个numa节点上都
void ** index_pemm_alloc_log(size_t size);

} // end namespace nnbtree