project(nnbtree)
cmake_minimum_required(VERSION 3.0)

# -march=native used to detect CPU flags and use corrosponding flush instruction
set(CMAKE_CXX_FLAGS "-march=native -Wall -O0 -g")

# Require C++17
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_STANDARD 17)

# nnbtree shared library
include_directories(
  include
)

set(NNBTREE_HEAD_LIST 
  include/
  include/nnbtree
)

set(NNBTREE_SRC
  src/index_btree.cc
)

add_library(nnbtree SHARED ${NNBTREE_SRC})
target_link_libraries(nnbtree pmem pmemobj pthread)
target_include_directories(nnbtree PUBLIC ${NNBTREE_HEAD_LIST})

# simple test
add_executable(simple_test tests/test.cc)
target_link_libraries(simple_test nnbtree rt m)

add_executable(simple_test_mixed tests/test.cc)
target_link_libraries(simple_test_mixed nnbtree rt m)
