#include <assert.h>
#include <error.h>
#include <limits.h>
#include <pthread.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define INFINITY 999999
#define SYSEXPECT(expr) do { if(!(expr)) { perror(__func__); exit(1); } } while(0)
#define error_exit(fmt, ...) do { fprintf(stderr, "%s error: " fmt, __func__, ##__VA_ARGS__); exit(1); } while(0);

int NVERTICES = -1;
int *GRAPH = NULL; // adjacency matrix
int *VERTICES = NULL;
int MIN_COST = INFINITY;
int *MIN_PATH = NULL;

// typedef struct node {
//     int cost;
//     int vertex;
//     int numVertices;
//     struct node *next;
// } node_t;

typedef struct tuple {
    int src;
    int des;
    int weight;
    struct tuple *next;
} tuple_t;

// node_t *graph_start; // adjacency list 
// node_t *graph_list;
// void init_adjacency_list(node_t **head, int cost, int vertex, int numVertices) {
//     *head = (node_t*)malloc(sizeof(node_t));
//     (*head)->cost = cost;
//     (*head)->vertex = vertex;
//     (*head)->numVertices = numVertices;
//     (*head)->next = NULL;
//     return;
// }

// void add_node(node_t **head, int cost, int vertex, int numVertices) {
//     node_t *new_node = (node_t*)malloc(sizeof(node_t));
//     if (*head == NULL) {
//         init_adjacency_list(head, cost, vertex, numVertices);
//         return;
//     } else {
//         new_node->next = NULL;
//         new_node->cost = cost;
//         new_node->vertex = vertex;
//         new_node->numVertices = numVertices;
//         new_node->next = *head;
//         *head = new_node;
//     }
//     return;
// }

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


// void print_graph(node_t **graph_list) {
//     for (int i = 0; i < NVERTICES; i++) {
//         node_t *temp_node = graph_list[i];
//         printf("%d", i);
//         while (temp_node) {
//             printf(" -- > %d ", temp_node->vertex);
//             temp_node = temp_node->next;
//         }
//         printf("\n");
//     }
// }

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

// inline static void set_edge(int i, int j, int value) {
//     int offset = i * NVERTICES + j;
//     VERTICES[offset] = value;
//     return;
// }

// inline static int get_edge(int i, int j) {
//     int offset = i * NVERTICES + j;
//     return VERTICES[offset];
// }

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
        for (int i = 0; i < NVERTICES; i++) {
            MIN_PATH[i] = dist[i];
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

    for (int i = 0; i < NVERTICES; i++) {
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

void bellman_matrix_start(int src, tuple_t *edges) {
    int* dist;
    dist = (int*)calloc(NVERTICES, sizeof(int));
    for (int i = 0; i < NVERTICES; i++) {
        dist[i] = INFINITY;
    }
    dist[src] = 0;
    // for (int i = 0; i < NVERTICES; i++) {
    //     for (int j = 0; j < NVERTICES; j++) {
    //         if (get_cost(i,j) != 0 && ((dist[i] + get_cost(i,j)) < dist[j])) {
    //             dist[j] = dist[i] + get_cost(i,j);
    //         }
    //     }
    // }

    for (int i = 0; i < NVERTICES-1; i++) {
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


void dijkstra_matrix_start(int src) {
    int dist[NVERTICES];

    bool sptSet[NVERTICES];

    for (int i = 0; i < NVERTICES; i++)
        dist[i] = INFINITY, sptSet[i] = false;

    // Distance of source vertex from itself is always 0
    dist[src] = 0;

    // Find shortest path for all vertices
    for (int count = 0; count < NVERTICES-1; count++) {
        int u = minDistance(dist, sptSet, NVERTICES);

    // Mark the picked vertex as processed
        sptSet[u] = true;

    // Update dist value of the adjacent vertices of the picked vertex.
        for (int v = 0; v < NVERTICES; v++) {
            if (!sptSet[v] && get_cost(u,v) && dist[u] != INFINITY 
                                        && dist[u]+get_cost(u,v) < dist[v])
            dist[v] = dist[u] + get_cost(u,v);
        }
    }

    // print the constructed distance array
    // printSolution(dist, NVERTICES);
    set_min(dist, NVERTICES);
}

void dijkstra_matrix_helper(int src, int* modifiedGraph) {
    int dist[NVERTICES];

    bool sptSet[NVERTICES];

    // Initialize all distances as INFINITE and stpSet[] as false
    for (int i = 0; i < NVERTICES; i++)
        dist[i] = INFINITY, sptSet[i] = false;

    dist[src] = 0;

    // Find shortest path for all vertices
    for (int count = 0; count < NVERTICES-1; count++) {
        int u = minDistance(dist, sptSet, NVERTICES);

        // Mark the picked vertex as processed
        sptSet[u] = true;

        // Update dist value of the adjacent vertices of the picked vertex.
        for (int vertex = 0; vertex < NVERTICES; vertex++) {
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

    for (int src = 0; src < NVERTICES; src++) {
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
    if (argc < 2) error_exit("Expecting one argument: [file name] \n");
    char *filename = argv[1];
    FILE *fp = fopen(filename, "r");
    if (fp == NULL) error_exit("Failed to open input file \"%s\"\n", filename);
    int scan_ret;
    scan_ret = fscanf(fp, "%d", &NVERTICES);
    if(scan_ret != 1) error_exit("Failed to read vertex count\n");
    GRAPH = (int*)calloc(NVERTICES * NVERTICES, sizeof(int));
    VERTICES = (int*)calloc(NVERTICES * NVERTICES, sizeof(int));
    MIN_PATH = (int*)calloc(NVERTICES, sizeof(int));
    SYSEXPECT(GRAPH != NULL);
    SYSEXPECT(VERTICES != NULL);
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

    // node_t *graph_list[NVERTICES];
    // for (int i = 0; i < NVERTICES; i++) {
    //     graph_list[i] = NULL;
    // }

    // for (int i = 0; i < NVERTICES; i++) {
    //     for (int j = 0; j < NVERTICES; j++) {
    //         if (get_cost(i, j) > 0) {
    //             add_node(&graph_list[i], get_cost(i, j), j, NVERTICES);
    //         }
    //     }
    // }
    
    // print_graph(graph_list);

    // uncomment this to run bellman-ford
    tuple_t *edges;
    edges = (tuple_t*)calloc((NVERTICES+1) * (NVERTICES+1), sizeof(tuple_t));
    for (int i = 0; i < NVERTICES; i++) {
        for (int j = 0; j < NVERTICES; j++) {
            if (get_cost(i,j) != 0) {
                add_tuple(&edges,i,j,get_cost(i,j));
            }
        }
    }

    struct timespec before, after;
    clock_gettime(CLOCK_REALTIME, &before);

    // run johnson 
    // johnson_matrix_start();

    // run dijkstra 
    // for (int src = 0; src < NVERTICES; src++) {
    //     dijkstra_matrix_start(src);
    // }

    // run bellman-ford
    for (int src = 0; src < NVERTICES; src++) {
        bellman_matrix_start(src, edges);
    }

    // clock 
    clock_gettime(CLOCK_REALTIME, &after);
    printSolution(MIN_PATH, NVERTICES);
    double delta_ms = (double)(after.tv_sec - before.tv_sec) * 1000.0 + (after.tv_nsec - before.tv_nsec) / 1000000.0;
    putchar('\n');
    printf("============ Time ============\n");
    printf("Time: %.3f ms (%.3f s)\n", delta_ms, delta_ms / 1000.0);
    printf("Min Cost: %d\n", MIN_COST);
    
    free(GRAPH);
    free(VERTICES);
    free(MIN_PATH);
    // wsp_print_result();
    return 0;
}
