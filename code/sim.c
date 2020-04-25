#include "sim.h"
#include "rutil.h"

#define MAX_MASS 10
#define MAXIMUM_RAD 100
#define KFACTOR 3
#define MFACTOR 3
#define MEPSILON 0.01



int choose_k (int rc) {
	return KFACTOR * rc * rc;
}

int choose_w (int M) {
	return (int)(MFACTOR * pow(M, 1 + MEPSILON));
}

bool equal_coord(coord_t a, coord_t b) {
  return a.x == b.x && a.y == b.y; 
}

cluster_t *init_graph(int radius) {
	int diameter = 2 * radius + 1;
	cluster_t *g = malloc(sizeof(cluster_t));
	if (g == NULL) return NULL;
	g->matrix = malloc(sizeof(char*) * diameter);
  g->radius = radius; 
  g->diameter = diameter; 
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
  float deg = next_random_float(seedp, 3.14);
  printf(" rb = %d || deg = %f\n", rb, deg);

  coord_t res;
  res.x = rb * sin(deg);
  res.y = rb * cos(deg);
  return res;
}

void create_walk(coord_t *walk, int rb, int k, coord_t start) {
  int i;
  random_t *seedp = malloc(sizeof(random_t));
  *seedp = DEFAULTSEED;
  for (i = 0; i < k; i++) {
    random_t f = next_random_float(seedp, 4);
    int dir = round(f + 1);
    // printf("dir = %d\n", dir);
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

void show_walk(coord_t *walk, int length) {
  int i ; 
  for (i = 0; i < length; i ++) {
    printf("walk i = %d : x = %d , y = %d \n", i, walk[i].x, walk[i].y); 
  }
}

coord_t **step_2 (param_t *params) {
	int i;
	int w = params->w;
	coord_t **all_walks = malloc(sizeof(coord_t*) * w);
	for (i = 0; i < w; i++){
		all_walks[i] = malloc(sizeof(coord_t) * params->k); // k steps
		create_walk(all_walks[i], params->rb, params->k, create_start(params->rb));
    show_walk(all_walks[i], params->k); 
	}
  return all_walks;
}

bool is_sticky(coord_t loc, cluster_t *cluster) {
  int locx = loc.x + cluster->radius; 
  int locy = loc.y + cluster->radius; 
  char **matrix = cluster->matrix; 
  return (matrix[locx - 1][locy] + matrix[locx + 1][locy] 
    + matrix[locx][locy-1] + matrix[locx][locy+1] > 0); 
}

int *step_3 (coord_t **walks, cluster_t *cluster, param_t *params) {
  int i, j;
  int *res = malloc(sizeof(int) * params->k);
  memset(res, -1, params->k * sizeof(int));
  for (i = 0; i < params->w; i++) {
    for (j = 0; j < params->k; j++) {
      coord_t loc = walks[i][j];
      if (is_sticky(loc, cluster)) {
        res[i] = j;
        break;
      }
    }
  }
  return res;
}

int step_4 (int* res, coord_t **walks, param_t *params) {
  int i, j, k;
  for (i = 0; i < params->w; i++) {
    for (j = 0; j < i; j++) {
      if (res[j] == -1)
        break;
      for (k = 0; k < params->k; k++) {
        if (equal_coord(walks[j][res[j]], walks[i][k]))
        {
          return i;
        }
          
      }
    }
  }
  return -1;
}

int step_5 (cluster_t *cluster, int* res, coord_t **walks, param_t *params, 
  int k, int *M) {
  int i;
  int rc = 0;
  k = k == -1 ? params->w : k;

  for (i = 0; i < k; i++) {
    if (res[i] == -1) {
      continue; 
    }

    coord_t tup = walks[i][res[i]];
    int this_rc = round(sqrt(tup.x * tup.x + tup.y * tup.y));
    if (cluster->matrix[tup.x + cluster->radius][tup.y + cluster->radius] == 0) {
      cluster->matrix[tup.x + cluster->radius][tup.y + cluster->radius] = 1;
      rc = this_rc > rc ? this_rc : rc;
      *M += 1; 
      printf("new particle = %d: x = %d, y = %d || rc = %d\n", i, tup.x, tup.y, rc);
    }
  }
  return cluster->radius < rc ? cluster->radius : rc;
}

void do_batch(cluster_t *cluster) {
	int rc = 1;
	int M = 1;
  int MAXITERS = 3; 
  int iters = 0;
  while (iters < MAXITERS) {
	// while (M < MAX_MASS) {
    printf("====== iter start (rc = %d) ======\n", rc); 
    param_t *params = step_1(rc, M);
    // printf("step 1 completed\n");
    coord_t **walks = step_2(params);
    // printf("step 2 completed\n");

    int *res = step_3(walks, cluster, params); 
    // printf("step 3 completed\n");

    int k = step_4(res,walks,params); 
    // printf("step 4 completed\n");

    int thisrc = step_5(cluster, res, walks, params, k, &M); 
    rc = thisrc > rc ? thisrc : rc ;
    printf("after rc = %d || M = %d\n", rc, M); 
    iters += 1;
	}
}



int main(int argc, char *argv[]) {
	cluster_t *g = init_graph(MAXIMUM_RAD);
	do_batch(g);
}

