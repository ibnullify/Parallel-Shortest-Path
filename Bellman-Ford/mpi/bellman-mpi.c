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
int *GRAPH = NULL; // adjacency matrix
int MIN_COST[2];
// MIN_COST[0] = INFINITY;
int *MIN_PATH = NULL;
int global_min[2];

int procID, nproc;

// tuple to store source node, dest node and edge weight
typedef struct tuple {
    int src;
    int des;
    int weight;
    struct tuple *next;
} tuple_t;


void init_tuple(tuple_t **head, int src, int des, int weight) {
    *head = (tuple_t*)malloc(sizeof(tuple_t));
    (*head)->src = src;
    (*head)->des = des;
    (*head)->weight = weight;
    (*head)->next = NULL;
    return;
}

void add_tuple(tuple_t **head, int src, int des, int weight) {
    tuple_t *new_tuple = (tuple_t*)malloc(sizeof(tuple_t));
    if (*head == NULL) {
        init_tuple(head, src, des, weight);
    }
    new_tuple->next = NULL;
    new_tuple->src = src;
    new_tuple->des = des;
    new_tuple->weight = weight;
    new_tuple->next = *head;
    *head = new_tuple;
}


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
        MIN_COST[0] = temp_sum;
        MIN_COST[1] = procID;
        int i;
        for (i = 0; i < NVERTICES+1; i++) {
            MIN_PATH[i] = dist[i];
        }
    }
    return;
}


void bellman_matrix_start(int src, tuple_t *edges) {
    int* dist;
    dist = (int*)calloc(NVERTICES, sizeof(int));
    int i;
    int chunk_size = NVERTICES / nproc;
    // initialize dist array each process initializes part
    int send[chunk_size];
    // for (i = 0; i < NVERTICES; i++) {
    //     dist[i] = INFINITY;
    // }
    for (i = 0; i < chunk_size; i++) {
        send[i] = INFINITY;
    }
    MPI_Allgather(send, chunk_size, MPI_INT, dist, chunk_size, MPI_INT, MPI_COMM_WORLD);
    dist[src] = 0;
    // int chunk_size = NVERTICES / nproc;
    // int send[chunk_size];
    for (i = 0; i < NVERTICES; i++) {
        // each process does a different edge 
        if (i % nproc != procID) continue;
        tuple_t *temp_head = edges;
        while (temp_head) { 
            if ((dist[temp_head->src] != INFINITY) && (dist[temp_head->src] + temp_head->weight < dist[temp_head->des])) {
                dist[temp_head->des] = dist[temp_head->src] + temp_head->weight;
            }
            temp_head = temp_head->next;
        }
    }
    set_min(dist, NVERTICES);
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

    tuple_t *edges;
    edges = (tuple_t*)calloc((NVERTICES+1) * (NVERTICES+1), sizeof(tuple_t));
    for (i = 0; i < NVERTICES; i++) {
        for (j = 0; j < NVERTICES; j++) {
            if (get_cost(i,j) != 0) {
                add_tuple(&edges,i,j,get_cost(i,j));
            }
        }
    }

    struct timespec before, after;
    clock_gettime(CLOCK_REALTIME, &before);


    // run bellman-ford
    // int src;
    // for (src = 0; src < NVERTICES-1; src++) {
    //     if (src % nproc != procID) continue;
    //     bellman_matrix_start(src, edges);
    // }
    bellman_matrix_start(0, edges);
    MPI_Allreduce(&MIN_COST, &global_min, 1, MPI_2INT, MPI_MINLOC, MPI_COMM_WORLD);
    // if (MIN_COST == global_min) printf("procID %d", procID);
    // clock 
    clock_gettime(CLOCK_REALTIME, &after);
    // 
    double delta_ms = (double)(after.tv_sec - before.tv_sec) * 1000.0 + (after.tv_nsec - before.tv_nsec) / 1000000.0;
    if (procID == global_min[1]) {
        printSolution(MIN_PATH, NVERTICES);
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
