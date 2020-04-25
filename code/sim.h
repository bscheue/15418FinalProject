#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>


typedef struct {

    int radius; // radius
    int diameter;

    bool **matrix; // matrix

} cluster_t;

typedef struct {
    int rb;
    int k;
    int w;
} param_t;

typedef struct {
    int i;
    int j;
} coord_t;

cluster_t *init_graph(int rc);

