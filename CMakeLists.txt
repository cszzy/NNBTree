project(nnbtree)
cmake_minimum_required(VERSION 3.0)

# -march=native used to detect CPU flags and use corrosponding flush instruction
set(CMAKE_CXX_FLAGS "-march=native -Wall")

# Require C++17
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CXX_STANDARD 17)

# Debug or Release
# Use command `cmake -DCMAKE_BUILD_TYPE=Release ..` to make a release build
if(CMAKE_BUILD_TYPE MATCHES Release)
  message(STATUS "Release build.")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
  add_definitions(-DNDEBUG)
else()
  message(STATUS "Debug build.")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
endif()

# nnbtree shared library
include_directories(
    ${PROJECT_SOURCE_DIR}/include/nnbtree
)

set(NNBTREE_SRC
    src/index_btree.cc
)

add_library(nnbtree SHARED ${NNBTREE_SRC})
target_link_libraries(nnbtree pmem pmemkv pthread)

# simple test
add_executable(simple_test tests/test.cc)
target_link_libraries(simple_test nnbtree)
