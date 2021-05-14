CXXFLAGS += -Wall -Wextra -pthread -fopenmp -g -O3 -DNDEBUG

.phony: all bellman-omp release

all: release

release: bellman-omp.c
	g++ bellman-omp.c -o bellman-omp $(CXXFLAGS)

clean:
	rm -f ./bellman-omp
	rm -f ./bellman-omp