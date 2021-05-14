CXXFLAGS += -Wall -Wextra -pthread -fopenmp -g -O3 -DNDEBUG

.phony: all johnson-mpi release

all: release

release: johnson-mpi.c
	mpicc johnson-mpi.c -o johnson-mpi $(CXXFLAGS)

clean:
	rm -f ./johnson-mpi
	rm -f ./johnson-mpi