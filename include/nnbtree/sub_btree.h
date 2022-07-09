#pragma once

// subtree分为在NVM和DRAM两种情况
    // 1. 初始时, 只有一个subtree, 即concurrent fast_fair
    // 2. 第一个subtree增长到一定高度, 则触发分裂. 
        // 此时上层树缓存到DRAM并从NVM删除(可以基于log, 但我觉得也可以阻塞请求，因为该操作只发生一次),
        // 底层作为subtree保留在NVM.
    // 3. 之后subtree可能仍然不断增长, 继续分裂(基于log)即可
    // 同时subtree也要进行缓存,淘汰和同步
    // 注意: fast_fair虽然实现了删除时的节点合并，但是为了性能实际并没有合并

#include "nnbtree.h"

namespace nnbtree {

// static pthread_mutex_t print_mtx;

/*
 * class SubTree
 */
SubTree::SubTree() {
  // root = (char *)new Page();
  // height = 1;
  sub_root = (char *)new(true) Page(PageType::NVM_SUBTREE_PAGE, 0); // 第一个subtree在pmem中
  clflush((char *)sub_root, sizeof(Page));
  height = 1;
  // clflush((char *)this, sizeof(SubTree));
  // printf("[nnbtree]: root is %p, SubTree is %p.\n", sub_root, this);
}

SubTree::SubTree(Page *root_) {
    if(root_ == nullptr) {
      sub_root = (char *)new(true) Page(PageType::NVM_SUBTREE_PAGE, 0);
      clflush((char *)sub_root, sizeof(Page));
      height = 1;
      // clflush((char *)this, sizeof(SubTree));
    } else {
      sub_root = (char *)root_;
      height = root_->GetLevel() + 1;
      // clflush((char *)this, sizeof(SubTree));
    }
    // printf("[nnbtree]: root is %p, SubTree is %p, height is %d.\n", sub_root, this, height);
}

void SubTree::setNewRoot(char *new_root) {
  this->sub_root = (char *)new_root;
  clflush((char *)&(this->sub_root), sizeof(char *));
  ++height;
  std::cout << "[subtree] setnewroot, height is " << height << std::endl;
}

char *SubTree::btree_search(entry_key_t key) {
  // std::lock_guard<std::mutex> l(subtree_lock);

  Page *p = (Page *)sub_root;

  while (p->hdr.leftmost_ptr != NULL) {
    Page *t = (Page *)p->linear_search(key);
    // p_assert(t != p->hdr.right_sibling_ptr, "should not happen");
    assert(t != p->hdr.right_sibling_ptr);
    if (t == p->hdr.right_sibling_ptr) {
      std::cout << "azheazhe" << std::endl;
    }
    p = t;
  }

  Page *t = NULL;
  while ((t = (Page *)p->linear_search(key)) == p->hdr.right_sibling_ptr) {
    // p_assert(false, "should not happen");
    // assert(false);
    std::cout << "azheazheazhe" << std::endl;
    p = t;
    if (!p) {
      break;
    }
  }

  if (!t) {
    printf("NOT FOUND %lu, t = %p\n", key, t);
    return NULL;
  }


  return (char *)t;
}

// insert the key in the leaf node
void SubTree::btree_insert(entry_key_t key, char *right) { // need to be string
  // std::lock_guard<std::mutex> l(subtree_lock);
  subtree_lock.lock();
retry:
  Page *p = (Page *)sub_root;

  while (p->hdr.leftmost_ptr != NULL) {
    Page *t = (Page *)p->linear_search(key);
    if (t == p->hdr.right_sibling_ptr) { // XXX : very important, 不能进入另一颗子树
      // p_assert(false, "should not happen");
      subtree_lock.unlock();
      return index_tree_root->btree_insert(key, right);
    }
    p = t;
  }

  if (!p->store(this, NULL, key, right, true, false)) { // store
    if (index_tree_root->has_hasindextree()) {
      subtree_lock.unlock();
      return index_tree_root->btree_insert(key, right);
    }
    // subtree_lock.unlock();
    // btree_insert(key, right);
    p_assert(false, "should not happen");
    goto retry;
  }

  subtree_lock.unlock();
}

// XXX: 需要判断store应该落在indextree还是subtree
// store the key into the node at the given level
void SubTree::btree_insert_internal(char *left, entry_key_t key, char *right,
                                  uint32_t level) {
  Page *p = (Page *)sub_root;

  if (level > (p->hdr.level))
    return;

  while (p->hdr.level > level)
    p = (Page *)p->linear_search(key);

  if (!p->store(this, NULL, key, right, true, false)) {
    btree_insert_internal(left, key, right, level);
  }
}

void SubTree::btree_delete(entry_key_t key) {
  // std::lock_guard<std::mutex> l(subtree_lock);
  subtree_lock.lock();
  
  Page *p = (Page *)sub_root;

  while (p->hdr.leftmost_ptr != NULL) { // 遍历到叶
    Page *t = (Page *)p->linear_search(key);
    p_assert(t != p->hdr.right_sibling_ptr, "should not happen");
    p = t;
  }

  Page *t = NULL;
  while ((t = (Page *)p->linear_search(key)) == p->hdr.right_sibling_ptr) {
    p_assert(false, "should not happen");
    p = t;
    if (!p)
      break;
  }

  if (p && t) {
    if (!p->remove(this, key, false, false)) {
      btree_delete(key);
    }
  } else {
    // printf("not found the key to delete %lu\n", key);
  }

  subtree_lock.unlock();
}

// 内部节点如果有目标key,则需要删除
void SubTree::btree_delete_internal(entry_key_t key, char *ptr, uint32_t level,
                                  entry_key_t *deleted_key,
                                  bool *is_leftmost_node, Page **left_sibling) {
  Page *p = (Page *)sub_root;

  if (level > (p->hdr.level))
    return;

  while (p->hdr.level > level) {
    p = (Page *)p->linear_search(key);
  }

  p->hdr.mtx->lock();

  if ((char *)p->hdr.leftmost_ptr == ptr) {
    *is_leftmost_node = true;
    p->hdr.mtx->unlock();
    return;
  }

  *is_leftmost_node = false;

  for (int i = 0; p->records[i].ptr != NULL; ++i) {
    if (p->records[i].ptr == ptr) {
      if (i == 0) {
        if ((char *)p->hdr.leftmost_ptr != p->records[i].ptr) { // 如果==,则根据FAST算法,该元素是无效的,忽略
          *deleted_key = p->records[i].key;
          *left_sibling = p->hdr.leftmost_ptr;
          p->remove(this, *deleted_key, false, false);
          break;
        }
      } else {
        if (p->records[i - 1].ptr != p->records[i].ptr) {
          *deleted_key = p->records[i].key;
          *left_sibling = (Page *)p->records[i - 1].ptr;
          p->remove(this, *deleted_key, false, false);
          break;
        }
      }
    }
  }

  p->hdr.mtx->unlock();
  return;
}

// Function to search keys from "min" to "max"
void SubTree::btree_search_range(entry_key_t min, entry_key_t max,
                               unsigned long *buf) {
  Page *p = (Page *)sub_root;

  while (p) {
    if (p->hdr.leftmost_ptr != NULL) {
      // The current Page is internal
      p = (Page *)p->linear_search(min);
    } else {
      // Found a leaf
      p->linear_search_range(min, max, buf);

      break;
    }
  }
}


void SubTree::btree_search_range(entry_key_t min, entry_key_t max, 
    std::vector<std::pair<uint64_t, uint64_t>> &result, int &size) {
    Page *p = (Page *)sub_root;

    while(p) {
        if(p->hdr.leftmost_ptr != NULL) {
        // The current Page is internal
            p = (Page *)p->linear_search(min);
        }
        else {
        // Found a leaf
            p->linear_search_range(min, max, result, size);
        break;
        }
    }
}

void SubTree::btree_search_range(entry_key_t min, entry_key_t max, void **values, int &size) {
    Page *p = (Page *)sub_root;

    while(p) {
        if(p->hdr.leftmost_ptr != NULL) {
        // The current Page is internal
            p = (Page *)p->linear_search(min);
        }
        else {
        // Found a leaf
            p->linear_search_range(min, max, values, size);
        break;
        }
    }
}

void SubTree::PrintInfo() {
    // printf("This is fast_fair b+ tree.\n");
    // printf("Node size is %lu, M path is %d.\n", sizeof(Page), cardinality);
    std::cout << "Tree height is " << height << std::endl << std::flush;
}

void SubTree::printAll() {
  // pthread_mutex_lock(&print_mtx);
  int total_keys = 0;
  Page *leftmost = (Page *)sub_root;
  printf("root: %p\n", leftmost);
  do {
    Page *sibling = leftmost;
    while (sibling) {
      if (sibling->hdr.level == 0) {
        total_keys += sibling->hdr.last_index + 1;
      }
      sibling->print();
      sibling = sibling->hdr.right_sibling_ptr;
    }
    printf("-----------------------------------------\n");
    leftmost = leftmost->hdr.leftmost_ptr;
  } while (leftmost);

  printf("total number of keys: %d\n", total_keys);
  // pthread_mutex_unlock(&print_mtx);
}

Page *SubTree::getRoot() {
  return (Page *)sub_root;
}

void SubTree::lock_subtree() {
  subtree_lock.lock();
}

void SubTree::unlock_subtree() {
  subtree_lock.unlock();
}

} // end namespace nnbtree