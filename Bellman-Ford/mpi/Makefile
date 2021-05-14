CXXFLAGS += -Wall -Wextra -pthread -fopenmp -g -O3 -DNDEBUG
LINKER += -I$(shell pwd)/mpi4.1/include -pthread -Wl,-rpath -Wl,$(shell pwd)/mpi4.1/lib -Wl,--enable-new-dtags -L$(shell pwd)/mpi4.1/lib -lmpi
.phony: all bellman-mpi release

all: release

release: bellman-mpi.c
	mpicc $(LINKER) bellman-mpi.c -o bellman-mpi $(CXXFLAGS)

clean:
	rm -f ./bellman-mpi
	rm -f ./bellman-mpi