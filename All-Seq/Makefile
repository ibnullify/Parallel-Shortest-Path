CXXFLAGS += -Wall -Wextra -pthread -fopenmp -g -O3 -DNDEBUG

.phony: all all-seq release

all: release

release: all-seq.c
	g++ all-seq.c -o all-seq $(CXXFLAGS)

clean:
	rm -f ./all-seq
	rm -f ./all-seq