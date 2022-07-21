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

#include <queue>

namespace nnbtree {

// static pthread_mutex_t print_mtx;

/*
 * class SubTree
 */
SubTree::SubTree(SubTreeStatus status) : left_sibling_subtree_(NULL), right_sibling_subtree_(NULL),
                       treelog_(NULL), tmp_hotness(0) {
  // root = (char *)new Page();
  // height_ = 1;
  sub_root_ = (char *)new(true) Page(PageType::NVM_SUBTREE_PAGE, 0); // 第一个subtree在pmem中
  clflush((char *)sub_root_, sizeof(Page));
  height_ = 1;
  
  nvm_root_ = sub_root_;
  subtree_status_ = status;
  memset(read_times_, 0, sizeof(read_times_));
  memset(write_times_, 0, sizeof(write_times_)); 
  memset(hotness_, 0, sizeof(hotness_));
  // clflush((char *)this, sizeof(SubTree));
  // printf("[nnbtree]: root is %p, SubTree is %p.\n", sub_root_, this);
}

SubTree::SubTree(Page *root_, SubTreeStatus status) : left_sibling_subtree_(NULL), right_sibling_subtree_(NULL),
                       treelog_(NULL), tmp_hotness(0) {
    if(root_ == nullptr) {
      sub_root_ = (char *)new(true) Page(PageType::NVM_SUBTREE_PAGE, 0);
      clflush((char *)sub_root_, sizeof(Page));
      height_ = 1;

      nvm_root_ = sub_root_;
      subtree_status_ = status;
      memset(read_times_, 0, sizeof(read_times_));
      memset(write_times_, 0, sizeof(write_times_));
      memset(hotness_, 0, sizeof(hotness_));
      // clflush((char *)this, sizeof(SubTree));
    } else {
      sub_root_ = (char *)root_;
      height_ = root_->GetLevel() + 1;

      nvm_root_ = sub_root_;
      subtree_status_ = status;
      memset(read_times_, 0, sizeof(read_times_));
      memset(write_times_, 0, sizeof(write_times_));
      memset(hotness_, 0, sizeof(hotness_));
      // clflush((char *)this, sizeof(SubTree));
    }
    // printf("[nnbtree]: root is %p, SubTree is %p, height_ is %d.\n", sub_root_, this, height_);
}

void SubTree::setNewRoot(char *new_root) {
  this->sub_root_ = (char *)new_root;
  clflush((char *)&(this->sub_root_), sizeof(char *));
  ++height_;
  std::cout << "[subtree] setnewroot, height_ is " << height_ << std::endl;
}

char *SubTree::btree_search(entry_key_t key) {
  // std::lock_guard<std::mutex> l(subtree_lock);

  Page *p = (Page *)sub_root_;

  while (p->hdr.leftmost_ptr != NULL) {
    // Page *t = (Page *)p->linear_search(key);
    // // p_assert(t != p->hdr.right_sibling_ptr, "should not happen");
    // assert(t != p->hdr.right_sibling_ptr);
    // if (t == p->hdr.right_sibling_ptr) {
    //   std::cout << "azheazhe" << std::endl;
    // }
    // p = t;

    p = (Page *)p->linear_search(key);
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

  read_times_[numa_map[my_thread_id]]++;

  if (!t) {
    printf("NOT FOUND %lu\n", key);
    return NULL;
  }


  return (char *)t;
}

// insert the key in the leaf node
void SubTree::btree_insert(entry_key_t key, char *right) { // need to be string
  // std::lock_guard<std::mutex> l(subtree_lock);
  subtree_lock.lock();
#ifdef CACHE_SUBTREE
  if (subtree_status_ == SubTreeStatus::IN_DRAM) {
    // treelog_->write_log(TreeLogType::INSERT, key, (uint64_t)right);
  } else if (subtree_status_ == SubTreeStatus::NEED_MOVE_TO_NVM) {
    // 刷回nvm
    move_to_nvm();
    // 释放日志，设置日志指针为null
    delete treelog_;
    treelog_ = nullptr;
  }
#endif
retry:
  Page *p = (Page *)sub_root_;

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
    p_assert(false, "should not happen");
    // subtree_lock.unlock();
    // btree_insert(key, right);
    
    goto retry;
  }

  write_times_[numa_map[my_thread_id]]++;

#ifdef CACHE_SUBTREE
  // 检测是否需要缓存
  if (subtree_status_ == SubTreeStatus::NEED_MOVE_TO_DRAM) {
    // 分配日志
    treelog_ = new TreeLog();
    assert(treelog_);

    // 移到内存
    move_to_dram();
  }
#endif

  subtree_lock.unlock();

  // if (index_tree_root->has_hasindextree());
  //   index_tree_root->btree_search(key);
}

// XXX: 需要判断store应该落在indextree还是subtree
// store the key into the node at the given level
void SubTree::btree_insert_internal(char *left, entry_key_t key, char *right,
                                  uint32_t level) {
  Page *p = (Page *)sub_root_;

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
  
  Page *p = (Page *)sub_root_;

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
  Page *p = (Page *)sub_root_;

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
  Page *p = (Page *)sub_root_;

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
    Page *p = (Page *)sub_root_;

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
    Page *p = (Page *)sub_root_;

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
    std::cout << "Tree height_ is " << height_ << std::endl << std::flush;
}

void SubTree::printAll() {
  // pthread_mutex_lock(&print_mtx);
  int total_keys = 0;
  Page *leftmost = (Page *)sub_root_;
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
  return (Page *)sub_root_;
}

void SubTree::lock_subtree() {
  subtree_lock.lock();
}

void SubTree::unlock_subtree() {
  subtree_lock.unlock();
}

void SubTree::move_to_nvm() {
  // 层序遍历dram tree, 不检测脏页了，直接写回算了
  Page *dram_root = (Page *)sub_root_;
  const Page *p = dram_root;
  p_assert(p, "[move to nvm]: root is null");
  Page *cur_page = new(true) Page(*dram_root);
  Page *new_root = cur_page;
  cur_page->hdr.page_type = PageType::NVM_SUBTREE_PAGE;

  std::queue<std::pair<const Page*, Page*>> page_queue;
  page_queue.push({p, cur_page});
  while (!page_queue.empty()) {
    int page_nums = page_queue.size();
    for (int i = 0; i < page_nums; i++) {
      auto page_pair = page_queue.front();
      page_queue.pop();
      p = page_pair.first;
      cur_page = page_pair.second;

      // copy child page
      Page *q = new(true) Page(*(p->hdr.leftmost_ptr));
      q->hdr.page_type = PageType::NVM_SUBTREE_PAGE;
      cur_page->hdr.leftmost_ptr = q;
      if (q->hdr.level != 0) // 非叶子节点
        page_queue.push({p->hdr.leftmost_ptr, q});

      for (int j = 0; j <= p->hdr.last_index; j++) {
        Page *t = new(true) Page(*(Page*)p->records[j].ptr);
        t->hdr.page_type = PageType::NVM_SUBTREE_PAGE;
        cur_page->records[j].ptr = (char *)t;
        if (t->hdr.level != 0) // 非叶子节点
          page_queue.push({(Page*)p->records[j].ptr, t});
        q->hdr.right_sibling_ptr = t;
        q = t;
      }
    }
  }

  // 修改子树根指针
  sub_root_ = (char *)new_root;
  
  // 获得前继子树锁，修改前继子树page指针
  if (left_sibling_subtree_ == nullptr)
    return;
retry:
  left_sibling_subtree_->lock_subtree();
  if (left_sibling_subtree_->right_sibling_subtree_ == this) { // 前继子树没有处于分裂的中间状态
    Page *p = (Page *)left_sibling_subtree_->sub_root_;
    Page *q = (Page *)this->sub_root_;
    for (int i = 0; i < MAX_SUBTREE_HEIGHT; i++) {
      p->hdr.right_sibling_ptr = q;
      p = (Page *)(p->records[p->hdr.last_index].ptr);
      q = q->hdr.leftmost_ptr;
    }
    assert(q == nullptr);
    left_sibling_subtree_->unlock_subtree();
    return;
  } else {
    SubTree *oldp = left_sibling_subtree_;
    left_sibling_subtree_ = left_sibling_subtree_->right_sibling_subtree_;
    oldp->unlock_subtree();
    goto retry;
  }

  // 设置nvm子树为新子树，释放旧子树
  Page *old_root = (Page *)nvm_root_;
  nvm_root_ = sub_root_;
  subtree_status_ = SubTreeStatus::IN_NVM;

  // 回收NVM空间
  std::queue<Page*> old_page_queue;
  old_page_queue.push(old_root);
  while (!old_page_queue.empty()) {
    int page_nums = old_page_queue.size();
    for (int i = 0; i < page_nums; i++) {
      p = old_page_queue.front();
      old_page_queue.pop();

      bool has_child = p->hdr.level != 0;
      if (has_child) // 还有孩子
        old_page_queue.push(p->hdr.leftmost_ptr);

      for (int j = 0; j <= p->hdr.last_index; j++) {
        if (has_child) // 还有孩子
          old_page_queue.push((Page*)p->records[j].ptr);
      }

      delete p; // 或者交给后台线程回收？
    }
  }

  // 回收DRAM空间
  std::queue<Page*> dram_page_queue;
  dram_page_queue.push(dram_root);
  while (!dram_page_queue.empty()) {
    int page_nums = dram_page_queue.size();
    for (int i = 0; i < page_nums; i++) {
      p = dram_page_queue.front();
      dram_page_queue.pop();

      // copy child page
      bool has_child = p->hdr.level != 0;
      if (has_child) // 还有孩子
        dram_page_queue.push(p->hdr.leftmost_ptr);

      for (int j = 0; j <= p->hdr.last_index; j++) {
        if (has_child) // 还有孩子
          dram_page_queue.push((Page*)p->records[j].ptr);
      }

      delete p; // XXX：如果有前台线程在对内存子树进行搜索，则会崩溃，需要延迟回收。
    }
  }
}

void SubTree::move_to_dram() {
  // 层序遍历subtree，copy到内存
  nvm_root_ = sub_root_;
  Page *p = (Page *)nvm_root_;
  p_assert(p, "[move to dram]: root is null");
  Page *dram_root = new(false) Page(*(Page *)nvm_root_); // 新生成的subtree的根page
  dram_root->hdr.page_type = PageType::DRAM_CACHETREE_PAGE;
  Page *cur_page = dram_root; // 当前操作的dram page
  std::queue<std::pair<Page *, Page*>> page_queue;
  page_queue.push({p, cur_page});
  while (!page_queue.empty()) {
    int page_nums = page_queue.size();
    for (int i = 0; i < page_nums; i++) {
      auto page_pair = page_queue.front();
      page_queue.pop();
      p = page_pair.first;
      cur_page = page_pair.second;

      // copy child page
      Page *q = new(false) Page(*(p->hdr.leftmost_ptr));
      q->hdr.page_type = PageType::DRAM_CACHETREE_PAGE;
      cur_page->hdr.leftmost_ptr = q;
      if (q->hdr.level != 0) // 还有孩子
        page_queue.push({p->hdr.leftmost_ptr, q});

      for (int j = 0; j <= p->hdr.last_index; j++) {
        Page *t = new(false) Page(*(Page*)p->records[j].ptr);
        t->hdr.page_type = PageType::DRAM_CACHETREE_PAGE;
        cur_page->records[j].ptr = (char *)t;
        if (t->hdr.level != 0) // 还有孩子
          page_queue.push({(Page*)p->records[j].ptr, t});
        q->hdr.right_sibling_ptr = t;
        q = t;
      }
    }
  }

  // 设置子树根，设置子树状态
  sub_root_ = (char*)dram_root;
  subtree_status_ = SubTreeStatus::IN_DRAM;

  // 获得前继subtree的锁（如果是第一颗子树则不需要获取了）
  // 修改前继subtree的page的right_sibling_ptr指向内存page，然后释放锁
  if (left_sibling_subtree_ == nullptr)
    return;
retry:
  left_sibling_subtree_->lock_subtree();
  if (left_sibling_subtree_->right_sibling_subtree_ == this) {
    Page *p = (Page *)left_sibling_subtree_->sub_root_;
    Page *q = (Page *)this->sub_root_;
    for (int i = 0; i < MAX_SUBTREE_HEIGHT; i++) {
      p->hdr.right_sibling_ptr = q;
      p = (Page *)(p->records[p->hdr.last_index].ptr);
      q = q->hdr.leftmost_ptr;
    }
    assert(q == nullptr);
    left_sibling_subtree_->unlock_subtree();
    return;
  } else {
    SubTree *oldp = left_sibling_subtree_;
    left_sibling_subtree_ = left_sibling_subtree_->right_sibling_subtree_;
    oldp->unlock_subtree();
    goto retry;
  }
}

void SubTree::cal_hotness() {
  // 当前不考虑numa，直接返回所有numa热度和作为总热度
  tmp_hotness = 0;
  // 现在热度计算就是简单地把各节点的读和写加在一起, 因为还不考虑numa
  for (int i = 0; i < numa_node_num; i++) {
    int cur_hotness = read_times_[i] + write_times_[i];
    read_times_[i] = write_times_[i] = 0;
    hotness_[i] = (cur_hotness + hotness_[i] >> 2) >> 1; // 瞎写的参数
    tmp_hotness += hotness_[i];
  }
}

uint64_t SubTree::get_hotness() const {
  return tmp_hotness;
}

void SubTree::setSubTreeStatus(SubTreeStatus status) {
  subtree_status_ = status;
}

SubTreeStatus SubTree::getSubTreeStatus() const {
  return subtree_status_;
}

} // end namespace nnbtree