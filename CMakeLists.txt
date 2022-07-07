project(nnbtree)
cmake_minimum_required(VERSION 3.0)

# -march=native used to detect CPU flags and use corrosponding flush instruction
set(CMAKE_CXX_FLAGS "-march=native -O3")

# Require C++17
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_STANDARD 17)

# nnbtree shared library
include_directories(
  include
  include/nnbtree
  include/fast_fair
)

set(NNBTREE_HEAD_LIST 
  include
  include/nnbtree
)

set(FAST_FAIR_HEAD_LIST 
  include
  include/fast_fair
)

set(NNBTREE_SRC
  src/index_btree.cc
  src/environment.cc
)

add_library(nnbtree SHARED ${NNBTREE_SRC})
target_link_libraries(nnbtree pmem pmemobj pthread)
target_include_directories(nnbtree PUBLIC ${NNBTREE_HEAD_LIST})

# simple test
add_executable(simple_test tests/test.cc)
target_link_libraries(simple_test nnbtree rt m)

# 如何增加编译选项-DMIXED?
add_executable(simple_test_mixed tests/test.cc)
target_link_libraries(simple_test_mixed nnbtree rt m)

# microbench
add_executable(microbench tests/micro_bench.cc)
target_link_libraries(microbench nnbtree)
target_include_directories(microbench PUBLIC ${FAST_FAIR_HEAD_LIST})
add_test(microbench microbench)

# multibench
add_executable(multibench tests/multi_bench.cc)
target_link_libraries(multibench nnbtree)
target_include_directories(multibench PUBLIC ${NNBTREE_HEAD_LIST})
add_test(multibench multibench)

#spinlock_test
add_executable(spinlock_test tests/spinlock_test.cc)
target_link_libraries(spinlock_test pthread)
add_test(spinlock_test spinlock_test)