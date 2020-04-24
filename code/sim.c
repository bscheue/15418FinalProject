#include "sim.h"
#include "rutil.h"

#define MAXIMUM_RAD 10
#define KFACTOR 3
#define MFACTOR 3
#define MEPSILON 0.01


int choose_k (int rc) {
	return KFACTOR * rc * rc;
}

int choose_w (int M) {
	return (int)(MFACTOR * pow(M, 1 + MEPSILON));
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

coord_t create_start(int rb) {
  random_t *seedp = malloc(sizeof(random_t));
  *seedp = DEFAULTSEED;
  float deg = next_random_float(seedp, 360);
  coord_t res;
  res.x = sin(deg);
  res.y = cos(deg);
  return res;
}

void create_walk(coord_t *walk, int rb, int k, coord_t start) {
  int i;
  random_t *seedp = malloc(sizeof(random_t));
  *seedp = DEFAULTSEED;
  for (i = 0; i < k; i++) {
    random_t f = next_random_float(seedp, 5);
    int dir = round(f);
    switch (dir) {
      case 1:
        walk[i].x = 1;
        walk[i].y = 0;
        break;
      case 2:
        walk[i].x = 0;
        walk[i].y = 1;
        break;
      case 3:
        walk[i].x = -1;
        walk[i].y = 0;
        break;
      case 4:
        walk[i].x = 0;
        walk[i].y = -1;
        break;
      default:
        printf("mistake\n");
    }
  }

  // compute prefix sum
  int accx = start.x;
  int accy = start.y;
  for (i = 0; i < k; i++) {
    accx += walk[i].x;
    accy += walk[i].y;
    walk[i].x = accx;
    walk[i].y = accy;
  }
}

param_t *step_1 (int rc, int M) {
	param_t *ret = malloc(sizeof(param_t));
	ret->rb = rc + 2;
	ret->k = choose_k(rc);
	ret->w = choose_w(M);
	return ret;
}

coord_t **step_2 (param_t *params) {
	int i;
	int w = params->w;
	coord_t **all_walks = malloc(sizeof(coord_t*) * w);
	for (i = 0; i < w; i++){
		all_walks[i] = malloc(sizeof(coord_t) * params->k); // k steps
		create_walk(all_walks[i], params->rb, params->k, create_start(params->rb));
	}
  return all_walks;
}

int *step_3 (coord_t **walks, cluster_t *cluster, param_t *params) {
  int i, j;
  int *res = malloc(sizeof(int) * params->k);
  memset(res, -1, params->k * sizeof(int));
  for (i = 0; i < params->w; i++) {
    for (j = 0; j < params->k; j++) {
      coord_t loc = walks[i][j];
      if (cluster->matrix[loc.x][loc.y]) {
        res[i] = j;
        break;
      }
    }
  }
  return res;
}

int step4 (int* res, coord_t **walks, param_t *params) {
  int i, j, k;
  for (i = 0; i < params->w; i++) {
    for (j = 0; j < i; j++) {
      if (res[j] == -1)
        break;
      for (k = 0; k < params->k; k++) {
        if (walks[j][res[j]].x == walks[i][k].x
            && walks[j][res[j]].y == walks[i][k].y)
          return i;
      }
    }
  }
  return -1;
}

int step5 (cluster_t *cluster, int* res, coord_t **walks, param_t *params, int k) {
  int i;
  int rc = 0;
  k = k == -1 ? params->w : k;

  for (i = 0; i < k; i++) {
    if (res[i] == -1)
      break;
    coord_t tup = walks[i][res[i]];
    int this_rc = round(sqrt(tup.x * tup.x + tup.y * tup.y));
    rc = this_rc > rc ? this_rc : rc;
    cluster->matrix[tup.x][tup.y] = 1;
  }
  return cluster->rc > rc ? cluster->rc : rc;
}

void do_batch(cluster_t *g) {
	int rc = 1;
	int M = 1;
	param_t *params = step_1(rc, M);
	step_2(params);
	// while (rc < GRAPH_RAD) {

	// }
}



int main(int argc, char *argv[]) {
	cluster_t *g = init_graph(MAXIMUM_RAD);
	do_batch(g);
}

