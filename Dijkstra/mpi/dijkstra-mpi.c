#include <assert.h>
#include <error.h>
#include <limits.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include "mpi.h"

#define INFINITY 999999
#define SYSEXPECT(expr) do { if(!(expr)) { perror(__func__); exit(1); } } while(0)
#define error_exit(fmt, ...) do { fprintf(stderr, "%s error: " fmt, __func__, ##__VA_ARGS__); exit(1); } while(0);

int NVERTICES = -1;
int NCORES = -1;
int *GRAPH = NULL; // adjacency matrix
// int *VERTICES = NULL;
int MIN_COST[2];
int *MIN_PATH = NULL;
int global_min[2];

int procID, nproc;

inline static void set_cost(int i, int j, int value) {
    assert(value > 0);
    int offset = i * NVERTICES + j;
    GRAPH[offset] = value;
    return;
}

inline static int get_cost(int i, int j) {
    int offset = i * NVERTICES + j;
    return GRAPH[offset];
}


inline static int get_value(int i, int j, int* matrix, int n) {
    int offset = i * n + j;
    return matrix[offset];
}

inline static void set_value(int i, int j, int value, int* matrix, int n) {
    int offset = i * n + j;
    matrix[offset] = value;
    return;
}

int minDistance(int dist[], int sptSet[], int numVertices) {
   // Initialize min value
    int min = INFINITY, min_index;
    min_index = 0;
    int v;
    for (v = 0; v < numVertices; v++)
        if (sptSet[v] == 0 && dist[v] <= min) {
            min = dist[v];
            min_index = v;
        }
   return min_index;
}

// A utility function to print the constructed distance array
int printSolution(int dist[], int n) {
    printf("Vertex   Distance from Source\n");
    int i;
    for (i = 0; i < n; i++)
        printf("%d \t\t %d\n", i, dist[i]);
    return 0;
}


int sum(int dist[], int n) {
    int sum = 0;
    int i;
    for (i = 0; i < n; i++) {
        if (dist[i] == INFINITY) {
            return INFINITY;
        }
        sum += dist[i];
    }
    return sum;
}

void set_min(int dist[], int n) {
    int temp_sum = sum(dist, n);
    if (temp_sum < MIN_COST[0]) {
        // MPI_Allreduce(&MIN_COST, &global_min, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);
        MIN_COST[0] = temp_sum;
        MIN_COST[1] = procID;
        int i;
        for (i = 0; i < NVERTICES; i++) {
            MIN_PATH[i] = dist[i];
        }
    }
    return;
}


void dijkstra_matrix_start(int src) {
    int dist[NVERTICES];
    int sptSet[NVERTICES];
    int i,count,v;
    for (i = 0; i < NVERTICES; i++)
        dist[i] = INFINITY, sptSet[i] = 0;

    dist[src] = 0;

    for (count = 0; count < NVERTICES-1; count++) {
        int u = minDistance(dist, sptSet, NVERTICES);
        sptSet[u] = 1;
        int chunk_size = NVERTICES / nproc;
        int send[chunk_size];
        // for (v = 0; v < NVERTICES; v++) {
        //     if (!sptSet[v] && get_cost(u,v) && dist[u] != INFINITY 
        //                                 && dist[u]+get_cost(u,v) < dist[v])
        //     dist[v] = dist[u] + get_cost(u,v);
        // }
        for (i = 0; i < chunk_size; i++) {
            // each process computes separate cost
            v = procID * chunk_size + i;
            if ((!sptSet[v]) && 
                (dist[v] > (dist[u] + get_cost(u, v))) && 
                get_cost(u,v)) {
                    send[i] = (dist[u] + get_cost(u, v));
                }
            else send[i] = dist[v];
        }
        // gather all computed costs 
        MPI_Allgather(send, chunk_size, MPI_INT, dist, chunk_size, MPI_INT, MPI_COMM_WORLD);
    }

    // print the constructed distance array
    // printSolution(dist, NVERTICES);
    set_min(dist, NVERTICES);
}



void dijkstra_list_start() {
    return;
}

int main(int argc, char **argv) {
    MPI_Init(NULL, NULL);
    MPI_Comm_rank(MPI_COMM_WORLD, &procID);
    MPI_Comm_size(MPI_COMM_WORLD, &nproc);
    if(argc < 2) error_exit("Expecting one argument: [file name]\n");
    char *filename = argv[1];
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) error_exit("Failed to open input file \"%s\"\n", filename);
    int scan_ret;
    scan_ret = fscanf(fp, "%d", &NVERTICES);
    if(scan_ret != 1) error_exit("Failed to read vertex count\n");
    GRAPH = (int*)calloc(NVERTICES * NVERTICES, sizeof(int));
    // VERTICES = (int*)calloc(NVERTICES * NVERTICES, sizeof(int));
    MIN_PATH = (int*)calloc(NVERTICES, sizeof(int));
    SYSEXPECT(GRAPH != NULL);
    // SYSEXPECT(VERTICES != NULL);
    SYSEXPECT(MIN_PATH != NULL);

    int i,j;
    for (i = 1; i < NVERTICES; i++) {
        for (j = 0; j < i; j++) {
            int t;
            scan_ret = fscanf(fp, "%d", &t);
            if (scan_ret != 1) error_exit("Failed to read cost(%d, %d)\n", i, j);
            set_cost(i, j, t); // comment out 
            set_cost(j, i, t);
        }
    }
    MIN_COST[0] = INFINITY;
    global_min[0] = INFINITY;
    struct timespec before, after;
    clock_gettime(CLOCK_REALTIME, &before);

    // run dijkstra 
    // int src;
    // for (src = 0; src < NVERTICES; src++) {
    //     if (src % nproc != procID) continue;
    //     dijkstra_matrix_start(src);
    // }
    dijkstra_matrix_start(0);
    MPI_Allreduce(&MIN_COST, &global_min, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);

    // clock 
    clock_gettime(CLOCK_REALTIME, &after);
    if (procID == global_min[1]) {
        printSolution(MIN_PATH, NVERTICES);
        double delta_ms = (double)(after.tv_sec - before.tv_sec) * 1000.0 + (after.tv_nsec - before.tv_nsec) / 1000000.0;
        putchar('\n');
        printf("============ Time ============\n");
        printf("Time: %.3f ms (%.3f s)\n", delta_ms, delta_ms / 1000.0);
        printf("Min Cost: %d\n", global_min[0]);
    }
    
    free(GRAPH);
    // free(VERTICES);
    free(MIN_PATH);
    // wsp_print_result();
    MPI_Finalize();
    return 0;
}
