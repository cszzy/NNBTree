#include "nvm_alloc.h"
#include "statistic.h"
#include "common_time.h"
#include "tree_log.h"
#include "numa_config.h"


namespace Common {
    std::map<std::string, Common::Statistic> timers;
    class Metcic g_metic;
    Stat stat;
}

namespace nnbtree {
    TreeLogPool *treelog_pool = nullptr;

    void init_treelogpool() {
        treelog_pool = new TreeLogPool();
        assert(treelog_pool);
        std::cout << "init treelog pool, pool nums:" << numa_node_num <<", log size: "
                 << LogSize / 1024 << "kB" << " per numa log nums: " << LogNum << std::endl;
    }
}

namespace NVM
{
Alloc *common_alloc = nullptr;
Alloc *data_alloc = nullptr;
Stat const_stat;
uint64_t  pmem_size = 0;

const size_t common_alloc_size = 100 * 1024 * 1024 * 1024UL;
const size_t data_alloc_size = 130 * 1024 * 1024 * 1024UL;

int env_init() {
    assert(!common_alloc);
    common_alloc = new  NVM::Alloc("/mnt/AEP0/zzy_common", common_alloc_size);
    assert(!data_alloc);
    data_alloc  = new  NVM::Alloc("/mnt/AEP0/zzy_data", data_alloc_size);
    return 0;
}

void env_exit()
{
    if(data_alloc) delete data_alloc;
    if(common_alloc) delete common_alloc;
    // if (nnbtree::treelog_pool) delete nnbtree::treelog_pool;
}

void show_stat()
{
    if(data_alloc)  data_alloc->Info();
    if(common_alloc)  common_alloc->Info();
}
} // namespace NVM
