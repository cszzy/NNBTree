// ycsb test

#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <future>
#include "ycsb/ycsb-c.h"

#include "nnbtree/sub_btree.h"
#include "fast_fair/fast_fair.h"
#include "nvm_alloc.h"

namespace nnbtree {

class fastfairDB : public ycsbc::KvDB
{
typedef uint64_t KEY_TYPE;
typedef uint64_t PAYLOAD_TYPE;
typedef FastFair::btree btree_t;

public:
fastfairDB() : tree_(nullptr) {}
fastfairDB(btree_t *btree) : tree_(btree) {}
virtual ~fastfairDB()
{
    delete tree_;
}

void Init()
{
    NVM::data_init();
    tree_ = new btree_t();
    NVM::pmem_size = 0;
}

void Info()
{
    std::cout << "NVM WRITE : " << NVM::pmem_size << std::endl;
    tree_->PrintInfo();
    NVM::show_stat();
}

void Bulk_load(const std::pair<uint64_t, uint64_t> data[], int size)
{
    for (int i = 0; i < size; i++)
    {
    tree_->btree_insert(data[i].first, (char *)data[i].second);
    }
}

void Close()
{
}
int Put(uint64_t key, uint64_t value)
{
    tree_->btree_insert(key, (char *)value);
    return 1;
}
int Get(uint64_t key, uint64_t &value)
{
    value = (uint64_t)tree_->btree_search(key);
    return 1;
}
int Update(uint64_t key, uint64_t value)
{
    tree_->btree_delete(key);
    tree_->btree_insert(key, (char *)value);
    return 1;
}

int Delete(uint64_t key)
{
    tree_->btree_delete(key);
    return 1;
}

int Scan(uint64_t start_key, int len, std::vector<std::pair<uint64_t, uint64_t>> &results)
{
    tree_->btree_search_range(start_key, UINT64_MAX, results, len);
    return 1;
}
void PrintStatic()
{
    NVM::show_stat();
}

private:
btree_t *tree_;
};

}