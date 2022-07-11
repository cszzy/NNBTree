#pragma once

#include "nnbtree.h"

// IndexTree
    // 和concurrent fast_fair一致，不同之处在于不需要clflush+mfence
namespace nnbtree {

// static pthread_mutex_t print_mtx;

/*
 * class IndexTree
 */

IndexTree::IndexTree(Page *page_, uint32_t level_) {
  root = (char *)page_;
  height = level_ + 1;
  std::cout << "[nnbtree]: indextree root is " << (void*)root << " , indextree is " << this << " , height is " << height << std::endl;
}


void IndexTree::setNewRoot(char *new_root) {
  this->root = new_root;
  ++height;
  std::cout << "[indextree] setnewroot, height is " << height << std::endl;
}

char *IndexTree::btree_search(entry_key_t key) {
  Page *p = (Page *)root;

  while (p->hdr.page_type != PageType::INDEXTREE_LAST_LEVEL_PAGE) {
    p = (Page *)p->linear_search(key);
  }

  SubTree *t = (SubTree*)(p->linear_search(key));
  while ((Page*)t == p->hdr.right_sibling_ptr) {
    // goto begin;
    p = (Page *)t;
    t = (SubTree*)(p->linear_search(key));
  }

  return t->btree_search(key);
}

// XXX: 搜索时indextree最后一层应该指向subtree，需要修改结构
// insert the key in the leaf node
void IndexTree::btree_insert(entry_key_t key, char *right) { // need to be string
// begin:
  Page *p = (Page *)root;

  while (p->hdr.page_type != PageType::INDEXTREE_LAST_LEVEL_PAGE) {
    p = (Page *)p->linear_search(key);
  }

  SubTree *t = (SubTree*)(p->linear_search(key));
  while ((Page*)t == p->hdr.right_sibling_ptr) {
    // goto begin;
    p = (Page *)t;
    t = (SubTree*)(p->linear_search(key));;
  }

  t->btree_insert(key, right);
}

// store the key into the node at the given level
void IndexTree::btree_insert_internal(char *left, entry_key_t key, char *right,
                                  uint32_t level) {
  Page *p = (Page *)root;
  if (level > (p->hdr.level))
    return;

  while (p->hdr.level > level)
    p = (Page *)p->linear_search(key);

  // p_assert(p->hdr.page_type == PageType::INDEXTREE_LAST_LEVEL_PAGE, "should not be false");
  if (!p->store(this, NULL, key, right, false, true)) {
    btree_insert_internal(left, key, right, level);
  }
}

void IndexTree::btree_delete(entry_key_t key) {
  Page *p = (Page *)root;

  while (p->hdr.page_type != PageType::INDEXTREE_LAST_LEVEL_PAGE) {
    p = (Page *)p->linear_search(key);
  }

  SubTree *t = (SubTree*)p->linear_search(key);

  t->btree_delete(key);
}

// 内部节点如果有目标key,则需要删除
void IndexTree::btree_delete_internal(entry_key_t key, char *ptr, uint32_t level,
                                  entry_key_t *deleted_key,
                                  bool *is_leftmost_node, Page **left_sibling) {
  Page *p = (Page *)root;
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
}

// Function to search keys from "min" to "max"
void IndexTree::btree_search_range(entry_key_t min, entry_key_t max,
                               unsigned long *buf) {
  Page *p = (Page *)root;

  while (p->hdr.page_type != PageType::INDEXTREE_LAST_LEVEL_PAGE) {
    p = (Page *)p->linear_search(min);
  }

  SubTree *t = (SubTree*)p->linear_search(min);

  t->btree_search_range(min, max, buf);
}


void IndexTree::btree_search_range(entry_key_t min, entry_key_t max, 
    std::vector<std::pair<uint64_t, uint64_t>> &result, int &size) {
    Page *p = (Page *)root;

    while (p->hdr.page_type != PageType::INDEXTREE_LAST_LEVEL_PAGE) {
      p = (Page *)p->linear_search(min);
    }

    SubTree *t = (SubTree*)p->linear_search(min);

    t->btree_search_range(min, max, result, size);
}

void IndexTree::btree_search_range(entry_key_t min, entry_key_t max, void **values, int &size) {
    Page *p = (Page *)root;

    while (p->hdr.page_type != PageType::INDEXTREE_LAST_LEVEL_PAGE) {
      p = (Page *)p->linear_search(min);
    }

    SubTree *t = (SubTree*)p->linear_search(min);

    t->btree_search_range(min, max, values, size);
}

Page * IndexTree::getRoot() {
    return (Page *)root;
}

void IndexTree::PrintInfo() {
    // printf("This is fast_fair b+ tree.\n");
    // printf("Node size is %lu, M path is %d.\n", sizeof(Page), cardinality);
    printf("Tree height is %d.\n", height);
}

void IndexTree::printAll() {
  // pthread_mutex_lock(&print_mtx);
  int total_keys = 0;
  Page *leftmost = (Page *)root;
  printf("root: %p\n", leftmost);
  do {
    if (leftmost->hdr.page_type == PageType::INDEXTREE_LAST_LEVEL_PAGE) {
      SubTree *t = (SubTree *)leftmost->records[0].ptr;
      leftmost = t->getRoot();
    }

    Page *sibling = leftmost;
    while (sibling) {
      if (sibling->hdr.level == 0) {
        total_keys += sibling->hdr.last_index + 1;
      }
      // sibling->print();
      sibling = sibling->hdr.right_sibling_ptr;
    }
    printf("-----------------------------------------\n");
    leftmost = leftmost->hdr.leftmost_ptr;
  } while (leftmost);

  printf("total number of keys: %d\n", total_keys);
  // pthread_mutex_unlock(&print_mtx);
}

} // end namespace nnbtree