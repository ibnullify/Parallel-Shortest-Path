CXXFLAGS += -Wall -Wextra -pthread -fopenmp -g -O3 -DNDEBUG
LINKER += -I$(shell pwd)/mpi4.1/include -pthread -Wl,-rpath -Wl,$(shell pwd)/mpi4.1/lib -Wl,--enable-new-dtags -L$(shell pwd)/mpi4.1/lib -lmpi
.phony: all dijkstra-mpi release

all: release

release: dijkstra-mpi.c
	mpicc dijkstra-mpi.c -o dijkstra-mpi $(CXXFLAGS)

clean:
	rm -f ./dijkstra-mpi
	rm -f ./dijkstra-mpi