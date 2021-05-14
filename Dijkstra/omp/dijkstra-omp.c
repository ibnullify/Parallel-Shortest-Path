#include <assert.h>
#include <error.h>
#include <limits.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <omp.h>

#define INFINITY 999999
#define SYSEXPECT(expr) do { if(!(expr)) { perror(__func__); exit(1); } } while(0)
#define error_exit(fmt, ...) do { fprintf(stderr, "%s error: " fmt, __func__, ##__VA_ARGS__); exit(1); } while(0);

int NVERTICES = -1;
int NCORES = -1;
int *GRAPH = NULL; // adjacency matrix
// int *VERTICES = NULL;
int MIN_COST = INFINITY;
int *MIN_PATH = NULL;


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

int minDistance(int dist[], bool sptSet[], int numVertices) {
   // Initialize min value
    int min = INFINITY, min_index;
    min_index = 0;
    for (int v = 0; v < numVertices; v++)
        if (sptSet[v] == false && dist[v] <= min) {
            min = dist[v];
            min_index = v;
        }
   return min_index;
}

// A utility function to print the constructed distance array
int printSolution(int dist[], int n) {
    printf("Vertex   Distance from Source\n");
    for (int i = 0; i < n; i++)
        printf("%d \t\t %d\n", i, dist[i]);
    return 0;
}


int sum(int dist[], int n) {
    int sum = 0;
    for (int i = 0; i < n; i++) {
        if (dist[i] == INFINITY) {
            return INFINITY;
        }
        sum += dist[i];
    }
    return sum;
}

void set_min(int dist[], int n) {
    int temp_sum = sum(dist, n);
    if (temp_sum < MIN_COST) {
        MIN_COST = temp_sum;
        // critical section for updating min_path
        #pragma omp critical 
        {
        for (int i = 0; i < NVERTICES; i++) {
            MIN_PATH[i] = dist[i];
        }
        }
    }
    return;
}


void dijkstra_matrix_start(int src) {
    int dist[NVERTICES];
    bool sptSet[NVERTICES];
    #pragma omp parallel for num_threads(NCORES)
    for (int i = 0; i < NVERTICES; i++)
        // paraellize initialization of distances 
        dist[i] = INFINITY, sptSet[i] = false;

    dist[src] = 0;

    for (int count = 0; count < NVERTICES-1; count++) {
        int u = minDistance(dist, sptSet, NVERTICES);
        sptSet[u] = true;
        #pragma omp parallel for num_threads(NCORES)
        for (int v = 0; v < NVERTICES; v++) {
            // parallelize over vertices 
            if (!sptSet[v] && get_cost(u,v) && dist[u] != INFINITY 
                                        && dist[u]+get_cost(u,v) < dist[v])
            dist[v] = dist[u] + get_cost(u,v);
        }
    }

    // print the constructed distance array
    // printSolution(dist, NVERTICES);
    set_min(dist, NVERTICES);
}



void dijkstra_list_start() {
    return;
}

int main(int argc, char **argv) {
    if(argc < 4 || strcmp(argv[1], "-p") != 0) error_exit("Expecting two arguments: -p [processor count] and [file name]\n");
    NCORES = atoi(argv[2]);
    if(NCORES < 1) error_exit("Illegal core count: %d\n", NCORES);
    char *filename = argv[3];
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


    for (int i = 1; i < NVERTICES; i++) {
        for (int j = 0; j < i; j++) {
            int t;
            scan_ret = fscanf(fp, "%d", &t);
            if (scan_ret != 1) error_exit("Failed to read cost(%d, %d)\n", i, j);
            set_cost(i, j, t); // comment out 
            set_cost(j, i, t);
        }
    }


    struct timespec before, after;
    clock_gettime(CLOCK_REALTIME, &before);

    // run dijkstra 
    // #pragma omp parallel for schedule(dynamic) num_threads(NCORES)
    dijkstra_matrix_start(0);
    // for (int src = 0; src < NVERTICES; src++) {
    //     dijkstra_matrix_start(src);
    // }


    // clock 
    clock_gettime(CLOCK_REALTIME, &after);
    printSolution(MIN_PATH, NVERTICES);
    double delta_ms = (double)(after.tv_sec - before.tv_sec) * 1000.0 + (after.tv_nsec - before.tv_nsec) / 1000000.0;
    putchar('\n');
    printf("============ Time ============\n");
    printf("Time: %.3f ms (%.3f s)\n", delta_ms, delta_ms / 1000.0);
    printf("Min Cost: %d\n", MIN_COST);
    
    free(GRAPH);
    // free(VERTICES);
    free(MIN_PATH);
    // wsp_print_result();
    return 0;
}
