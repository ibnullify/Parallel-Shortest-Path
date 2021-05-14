CXXFLAGS += -Wall -Wextra -pthread -fopenmp -g -O3 -DNDEBUG

.phony: all dijkstra-omp release

all: release

release: dijkstra-omp.c
	g++ dijkstra-omp.c -o dijkstra-omp $(CXXFLAGS)

clean:
	rm -f ./dijkstra-omp
	rm -f ./dijkstra-omp