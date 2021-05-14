CXXFLAGS += -Wall -Wextra -pthread -fopenmp -g -O3 -DNDEBUG

.phony: all wsp release

all: release

release: johnson-omp.c
	g++ johnson-omp.c -o johnson-omp $(CXXFLAGS)

clean:
	rm -f ./omp
	rm -f ./omp