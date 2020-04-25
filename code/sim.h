#include <stdlib.h>
#include <stdint.h>
#include <math.h>


typedef struct {

    int radius; // radius 
    int diameter; 

    char **matrix; // matrix

} cluster_t;

typedef struct {
    int rb; 
    int k; 
    int w; 
} param_t; 

typedef struct {
    int x;
    int y; 
} coord_t; 

cluster_t *init_graph(int rc); 

