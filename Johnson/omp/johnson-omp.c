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
int *GRAPH = NULL; // adjacency matrix
// int *VERTICES = NULL;
int MIN_COST = INFINITY;
int *MIN_PATH = NULL;
int NCORES = -1;


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
        // critical section update shared min_path 
        #pragma omp critical 
        {
        for (int i = 0; i < NVERTICES; i++) {
            MIN_PATH[i] = dist[i];
        }
        }
    }
    return;
}

int* bellman_matrix_helper(tuple_t* edges) {
    // initialize dist 
    int* dist;
    dist = (int*)calloc(NVERTICES+1, sizeof(int));
    for (int i = 0; i < NVERTICES+1; i++) {
        dist[i] = INFINITY;
    }
    for (int i = 0; i < NVERTICES; i++) {
        add_tuple(&edges, NVERTICES, i, 0);
        // set_value(NVERTICES, i, 0, edges, NVERTICES+1);
    }
    dist[NVERTICES] = 0;
    #pragma omp parallel for schedule(dynamic) num_threads(NCORES)
    for (int i = 0; i < NVERTICES; i++) {
        // each thread takes a separate vertex
        tuple_t *temp_head = edges;
        while (temp_head) {
            if ((temp_head->src != INFINITY) && (dist[temp_head->src] + temp_head->weight < dist[temp_head->des])) {
                dist[temp_head->des] = dist[temp_head->src] + temp_head->weight;
            }
            temp_head = temp_head->next;
        }
    }
    return dist;
}



void dijkstra_matrix_helper(int src, int* modifiedGraph) {
    int dist[NVERTICES];

    bool sptSet[NVERTICES];

    for (int i = 0; i < NVERTICES; i++)
        dist[i] = INFINITY, sptSet[i] = false;

    dist[src] = 0;

    for (int count = 0; count < NVERTICES-1; count++) {
        int u = minDistance(dist, sptSet, NVERTICES);

        // Mark the picked vertex as processed
        sptSet[u] = true;

        // Update dist value of the adjacent vertices of the picked vertex.
        #pragma omp parallel for schedule(dynamic) num_threads(NCORES)
        for (int vertex = 0; vertex < NVERTICES; vertex++) {
            // each thread takes a separate vertex
            if ((!sptSet[vertex]) && 
                (dist[vertex] > (dist[u] + get_value(u, vertex, modifiedGraph, NVERTICES))) && 
                get_cost(u,vertex)) {
                    dist[vertex] = (dist[u] + get_value(u, vertex, modifiedGraph, NVERTICES));
                }
        }
    }

    // print the constructed distance array
    // printSolution(dist, NVERTICES);
    set_min(dist, NVERTICES);
}

void johnson_matrix_start() {
    tuple_t *edges;
    edges = (tuple_t*)calloc((NVERTICES+1) * (NVERTICES+1), sizeof(tuple_t));
    for (int i = 0; i < NVERTICES; i++) {
        for (int j = 0; j < NVERTICES; j++) {
            if (get_cost(i,j) != 0) {
                add_tuple(&edges,i,j,get_cost(i,j));
            }
        }
    }

    int* modifyWeights;
    modifyWeights = (int*)calloc(NVERTICES+1, sizeof(int));
    modifyWeights = bellman_matrix_helper(edges);
    int* modifiedGraph = NULL;
    modifiedGraph = (int*)calloc(NVERTICES * NVERTICES, sizeof(int));
    for (int i = 0; i < NVERTICES; i++) {
        for (int j = 0; j < NVERTICES; j++) {
            set_value(i, j, 0, modifiedGraph, NVERTICES);
        }
    }

    for (int i = 0; i < NVERTICES; i++) {
        for (int j = 0; j < NVERTICES; j++) {
            if (get_cost(i,j) != 0) {
                int value = get_cost(i,j) + modifyWeights[i] - modifyWeights[j];
                set_value(i, j, value, modifiedGraph, NVERTICES);
            }
        }
    }
    #pragma omp parallel for schedule(dynamic) num_threads(NCORES)
    for (int src = 0; src < NVERTICES; src++) {
        // call dijkstras on different threads 
        dijkstra_matrix_helper(src, modifiedGraph);
    }
    free(modifiedGraph);
    free(modifyWeights);
    free(edges);
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

    // run johnson 
    johnson_matrix_start();


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
