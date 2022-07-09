#include "numa_config.h"
#include <iostream>

namespace nnbtree {

pmem::obj::pool_base pop_numa[nnbtree::numa_max_node];
int8_t numa_map[nnbtree::max_thread_num];
thread_local int my_thread_id;

void bindCore(uint16_t core) {
    // printf("bind to %d\n", core);
    cpu_set_t cpuset;
    CPU_ZERO(&cpuset);
    CPU_SET(core, &cpuset);
    int rc = pthread_setaffinity_np(pthread_self(), sizeof(cpu_set_t), &cpuset);
    if (rc != 0) {
        printf("can't bind core %d!", core);
        exit(-1);
    }
}

void init_numa_pool() {
    // core-thread_id  => numa_id
    for (int i = 0; i < 32; i++) {
        numa_map[i] = 0;
    }

    for (int i = 32; i < 64; i++) {
        numa_map[i] = 1;
    }

    for (int i = 0; i < nnbtree::numa_node_num; i++) {
        std::string pool_name = std::string("/mnt/AEP") 
                                + std::to_string(i) + "/numa";
        std::cout << "numa " << i << " pool: " << pool_name;
        remove(pool_name.c_str());
        pop_numa[i] = pmem::obj::pool<int>::create(
            pool_name, "WQ", PMEMOBJ_MIN_POOL * 8 * 1024, S_IWUSR | S_IRUSR
        );
    }
}

void *index_pmem_alloc(size_t size) {
    PMEMoid oid;
    if (pmemobj_alloc(pop_numa[numa_map[my_thread_id]].handle(), 
            &oid, size, 0, nullptr, nullptr)) {
        fprintf(stderr, "fail to alloc nvm\n");
        exit(-1);
    }

    return (void *)pmemobj_direct(oid);
}

void index_pmem_free(void *ptr)
{
  auto f_oid = pmemobj_oid(ptr);
  pmemobj_free(&f_oid);
}

}

