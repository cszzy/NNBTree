#pragma once

#include "nnbtree.h"
#include <unordered_set>
namespace nnbtree {

class Statistics {
    public:
        Statistics() {}

        ~Statistics() {}

        void insert_subtree(SubTree * const subtree) {
            all_subtree_lock.lock();
            all_subtree.push_back(subtree);
            all_subtree_lock.unlock();
        }

        // 选择前k个热点子树，目前只考虑缓存
        void select_topk(int k) {  
            // 把子树复制过来
            all_subtree_lock.lock();
            if (staticstic_subtree.size() == 0) {
                staticstic_subtree.assign(all_subtree.begin(), all_subtree.end());
            } else {
                staticstic_subtree.insert(staticstic_subtree.end(), 
                        all_subtree.begin() + staticstic_subtree.size(), all_subtree.end());
            }
            all_subtree_lock.unlock();


            // 计算所有子树热度
            for (auto subtree : staticstic_subtree) {
                subtree->cal_hotness();
            }

            // O(N)得到topk，同时将热点nvm子树设置为待缓存，将非热点nvm子树设置为待淘汰
            topk(staticstic_subtree, k);
        } 

        // 

    private:
        void topk(std::vector<SubTree *> &arr, int k) {
            if (arr.size() < k) {
                // subtree总数小于k，不进行内存缓存
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
                        } else if (arr[i]->getSubTreeStatus() == SubTreeStatus::NEED_MOVE_TO_NVM) {
                            arr[i]->setSubTreeStatus(SubTreeStatus::NEED_MOVE_TO_DRAM);
                            std::cout << "hasn't been move to nvm" << std::endl;
                        }
                    }

                    for (auto iter = topk_subtree.begin(); iter != topk_subtree.end(); iter++) {
                        (*iter)->setSubTreeStatus(SubTreeStatus::NEED_MOVE_TO_NVM);
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
            int pos = low + random() % (high - low + 1);
            // int pos = (low + high) >> 1;
            uint64_t pivot = arr[pos]->get_hotness();
            std::swap(arr[low], arr[pos]);
            int left = low + 1;
            int right = high;
            while (true) {
                while (left <= right && arr[left]->get_hotness() <= pivot) {
                    left++;
                }

                while (left <= right && arr[right]->get_hotness() >= pivot) {
                    right--;
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

extern Statistics *statis_;

} // end namespace nnbtree