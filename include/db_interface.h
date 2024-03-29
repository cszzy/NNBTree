// ycsb test

#pragma once

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <future>
#include "ycsb/ycsb-c.h"

#include "nnbtree/nnbtree.h"
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
    tree_ = new btree_t();
#ifdef TEST_PMEM_SIZE
    NVM::pmem_size = 0;
#endif
}

void Info()
{
#ifdef TEST_PMEM_SIZE
    std::cout << "NVM WRITE : " << NVM::pmem_size << std::endl;
#endif
    tree_->PrintInfo();
    NVM::show_stat();
    // uint64_t space = 0;
    // tree_->CalculateSapce(space);
    // std::cout << "calculate space: " << space * 1.0 / 1024 / 1024/ 1024 << " GB" << std::endl;
    // tree_->levelTraverse();
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

class nnbtreeDB : public ycsbc::KvDB
{
typedef uint64_t KEY_TYPE;
typedef uint64_t PAYLOAD_TYPE;
typedef nnbtree::NNBTree btree_t;

public:
nnbtreeDB() : tree_(nullptr) {}
nnbtreeDB(btree_t *btree) : tree_(btree) {}
virtual ~nnbtreeDB()
{
    delete tree_;
}

void Init()
{
    tree_ = new btree_t();
#ifdef TEST_PMEM_SIZE
    NVM::pmem_size = 0;
#endif
    init_treelogpool();
}

void Info()
{
#ifdef TEST_PMEM_SIZE
    std::cout << "NVM WRITE : " << NVM::pmem_size << std::endl;
#endif
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
    value = (uint64_t)(tree_->btree_search(key));
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