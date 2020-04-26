#include "sim.h"
#include "rutil.h"

#define MAX_MASS 10
#define MAXIMUM_RAD 30
#define KFACTOR 3
#define MFACTOR 3
#define MEPSILON 0.01

// maximum walk length
int choose_k(int rc) { return KFACTOR * rc * rc; }

// number of walkers
int choose_w(int M) { return (int)(MFACTOR * pow(M, 1 + MEPSILON)); }

bool equal_coord(coord_t a, coord_t b) { return a.i == b.i && a.j == b.j; }

cluster_t *init_graph(int radius) {
  int diameter = 2 * radius + 1;
  cluster_t *g = malloc(sizeof(cluster_t));
  if (g == NULL)
    return NULL;
  g->matrix = malloc(sizeof(bool*) * diameter);
  g->radius = radius;
  g->diameter = diameter;
  int i;
  for (i = 0; i < diameter; i++) {
    g->matrix[i] = calloc(diameter, sizeof(bool));
  }
  g->matrix[radius][radius] = 1;

  return g;
}

coord_t create_start(int rb, random_t *seed) {
  float deg = next_random_float(seed, 3.1415);
  coord_t res;
  printf("resi before = %d, j = %d\n", res.i, res.j);
  res.i = (int)(rb * sin(deg));
  res.j = (int)(rb * cos(deg));
  printf("=== resi = %d, resj = %d, resi calc = %d, resj calc = %d\n", res.i, res.j, (int)(rb * sin(deg)), (int)(rb * cos(deg)));
  return res;
}

void create_walk(coord_t *walk, int rb, int k, coord_t start, random_t *seed) {
  int i;
  for (i = 1; i < k; i++) {
    random_t f = next_random_float(seed, 4);
    int dir = round(f + 1);
    switch (dir) {
    case 1:
      walk[i].i = 1;
      walk[i].j = 0;
      break;
    case 2:
      walk[i].i = 0;
      walk[i].j = 1;
      break;
    case 3:
      walk[i].i = -1;
      walk[i].j = 0;
      break;
    case 4:
      walk[i].i = 0;
      walk[i].j = -1;
      break;
    default:
      printf("mistake\n");
    }
  }

  // compute prefix sum
  int acci = start.i;
  int accj = start.j;
  walk[0].i = acci;
  walk[0].j = accj;
  for (i = 1; i < k; i++) {
    acci += walk[i].i;
    accj += walk[i].j;
    walk[i].i = acci;
    walk[i].j = accj;
  }
}

param_t *step_1(int rc, int M) {
  param_t *ret = malloc(sizeof(param_t));
  ret->rb = rc + 2;
  ret->k = choose_k(rc);
  ret->w = choose_w(M);
  return ret;
}

void show_walk(int idx, coord_t *walk, int length) {
  /* int i; */
  /* printf("walk i = %d : ", idx); */
  /* for (i = 0; i < length; i++) { */
     /* printf("(%d, %d),", walk[i].i, walk[i].j); */
  /* } */
  /* printf("\n"); */
}

coord_t **step_2(param_t *params, random_t *seeds) {
  int i;
  int w = params->w;
  coord_t **all_walks = malloc(sizeof(coord_t *) * w);
  for (i = 0; i < w; i++) {
    all_walks[i] = malloc(sizeof(coord_t) * params->k); // k steps
    create_walk(all_walks[i], params->rb, params->k,
                create_start(params->rb, seeds + i), seeds + i);
    show_walk(i, all_walks[i], params->k);
  }
  return all_walks;
}

bool is_sticky(coord_t loc, cluster_t *cluster) {
  int locx = loc.i + cluster->radius;
  int locy = loc.j + cluster->radius;
  bool **matrix = cluster->matrix;
  /* printf("diam %d\n", cluster->diametr); */
  return
    locx >= 1 && locx < cluster->diameter - 1
    && locy >= 1 && locy < cluster->diameter - 1
    && (matrix[locx - 1][locy]
    ||  matrix[locx + 1][locy]
    ||  matrix[locx][locy - 1]
    ||  matrix[locx][locy + 1]);
}

int *step_3(coord_t **walks, cluster_t *cluster, param_t *params) {
  int i, j;
  int *res = malloc(sizeof(int) * params->w);
  memset(res, -1, params->w * sizeof(int));
  for (i = 0; i < params->w; i++) {
    for (j = 0; j < params->k; j++) {
      coord_t loc = walks[i][j];
      if (is_sticky(loc, cluster)) {
        res[i] = j;
        /* printf("^ particle %d is sticky at step %d \n", i, j); */
        break;
      }
    }
  }
  return res;
}

int step_4(int *res, coord_t **walks, param_t *params) {
  int i, j, k;
  for (i = 0; i < params->w; i++) {
    for (j = 0; j < i; j++) {
      if (res[j] == -1)
        break;
      for (k = 0; k < params->k; k++) {
        if (equal_coord(walks[j][res[j]], walks[i][k])) {
          /* printf("interference at particle %d\n", i); */
          return i;
        }
      }
    }
  }
  /* printf("NO INTERFERENCE\n"); */

  return -1;
}

int step_5(cluster_t *cluster, int *res, coord_t **walks, param_t *params,
           int k, int *M) {
  int i;
  int rc = 0;
  k = k == -1 ? params->w : k;

  for (i = 0; i < k; i++) {
    if (res[i] == -1) {
      continue;
    }

    coord_t tup = walks[i][res[i]];
    int this_rc = round(sqrt(tup.i * tup.i + tup.j * tup.j));
    if (cluster->matrix[tup.i + cluster->radius][tup.j + cluster->radius] == 0) {
      cluster->matrix[tup.i + cluster->radius][tup.j + cluster->radius] = true;
      rc = this_rc > rc ? this_rc : rc;
      *M += 1;
      printf("new particle = %d: i = %d, j = %d || rc = %d || resi = %d\n", i, tup.i, tup.j,
             rc, res[i]);
    }
  }
  return cluster->radius < rc ? cluster->radius : rc;
}

void view_cluster(cluster_t *g) {
  int i, j;
  for (i = 0; i < g->diameter; i++) {
    for (j = 0; j < g->diameter; j++) {
      if (i == g->radius && j == g->radius) printf("C ");
      else if (g->matrix[i][j] == 0) printf("- ");
      else printf("+ ");
    }
    printf("\n");
  }
}


void do_batch(cluster_t *cluster, random_t *seeds) {
  int rc = 1;
  int M = 1;
  int MAXITERS = 15;
  int iters = 0;
  while (iters < MAXITERS) {
    // while (M < MAX_MASS) {
    printf("====== iter start (rc = %d) ======\n", rc);
    param_t *params = step_1(rc, M);
    // printf("step 1 completed\n");
    coord_t **walks = step_2(params, seeds);
    // printf("step 2 completed\n");

    int *res = step_3(walks, cluster, params);
    // printf("step 3 completed\n");

    int k = step_4(res, walks, params);
    // printf("step 4 completed\n");

    int thisrc = step_5(cluster, res, walks, params, k, &M);
    rc = thisrc > rc ? thisrc : rc;
    /* printf("after rc = %d || M = %d\n", rc, M); */
    iters++;
    /* view_cluster(cluster); */
  }
}

random_t *init_seeds() {
  int seeds_size = 100000;
  random_t *seeds = malloc(sizeof(random_t) * seeds_size);
  random_t *seedseed = malloc(sizeof(random_t));
  *seedseed = DEFAULTSEED;
  int i;
  for (i = 0; i < seeds_size; i++) {
    seeds[i] = (random_t)next_random_float(seedseed, seeds_size);
  }
  return seeds;
}

int main(int argc, char *argv[]) {
  cluster_t *g = init_graph(MAXIMUM_RAD);
  random_t *seeds = init_seeds();
  do_batch(g, seeds);
  view_cluster(g);
  return 0;
}
