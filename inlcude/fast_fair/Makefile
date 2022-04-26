.PHONY: all clean
.DEFAULT_GOAL := all

LIBS=-lrt -lm -lpthread -lpmem -lpmemobj
CFLAGS=-O0 -std=c++17 -g 

output = btree_concurrent btree_concurrent_mixed

all: main

main: test.cpp
	g++ $(CFLAGS) -o btree_concurrent test.cpp $(LIBS) -DCONCURRENT
	g++ $(CFLAGS) -o btree_concurrent_mixed test.cpp $(LIBS) -DCONCURRENT -DMIXED

clean: 
	rm $(output)