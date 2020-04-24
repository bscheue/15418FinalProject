#include "sim.h"
#include "rutil.h"
#define MAXIMUM_RAD 10
#define KFACTOR 3
#define MFACTOR 3
#define MEPSILON 0.01

int main(int argc, char *argv[]) {
	cluster_t *g = init_graph(MAXIMUM_RAD); 
	do_batch(g); 
}

cluster_t *init_graph(int radius) {
	int diameter = 2 * radius + 1; 
	cluster_t *g = malloc(sizeof(cluster_t)); 
	if (g == NULL) return NULL; 
	g->matrix = malloc(sizeof(char*) * diameter); 
	int i; 
	for (i = 0; i < diameter; i ++) {
		g->matrix[i] = calloc(diameter, sizeof(char)); 
	}
	g->matrix[radius][radius] = 1; 

	return g; 
}

param_t *step_1 (int rc, int M) {
	param_t *ret = malloc(sizeof(param_t)) ; 
	ret->rb = rc + 2; 
	ret->k = choose_k(rc); 
	ret->w = choose_w(M); 
	return ret; 
}

coord_t **step_2 (param_t *params) {
	int i ; 
	int w = params->w; 
	coord_t **all_walks = malloc(sizeof(coord_t*) * w); 
	for (i = 0; i < w; i ++ ){
		all_walks[i] = malloc(sizeof(coord_t) * params->k); // k steps 
		create_walk(all_walks[i], params->rb, params->k, create_start(params->rb)); 
	}
}


void do_batch(cluster_t *g) {
	int rc = 1; 
	int M = 1; 
	param_t *params = step_1(rc, M); 
	step_2(params); 
	// while (rc < GRAPH_RAD) {

	// }
}

void create_walk(coord_t *walk, int rb, int k, coord_t* carray ) {
	int i = 0; 


}

coord_t create_start(int rb) {

}

int choose_k (int rc) {
	return KFACTOR * rc * rc; 
}

int choose_w (int M) {
	return (int)(MFACTOR * pow(M, 1 + MEPSILON)); 
}