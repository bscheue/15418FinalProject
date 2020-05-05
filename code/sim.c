#include <stdlib.h>
#include <stdint.h>
#include <math.h>
#include <stdbool.h>

#include "rutil.h"
#include "image.h"
#include "cycletimer.h"
#include "instrument.h"

#ifndef OMP
#define OMP 0
#endif

#if OMP
#include <omp.h>
#else
#include "fake_omp.h"
#endif

#define TAU 6.28
#define MAX_MASS 5000
#define MAX_RAD 300
#define RAD_THRESHOLD 75
#define KFACTOR 10
#define MFACTOR 0.1
#define MEPSILON 0.01
#define MIN_WALKERS 5
#define MAX_ITERS 300

#define TRACKING true


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

// maximum walk length
int choose_k(int rc) { return KFACTOR * rc * rc; }

// number of walkers
int choose_w(int M) { return MIN_WALKERS + MFACTOR * pow(M, 1 + MEPSILON); }

bool equal_coord(coord_t a, coord_t b) { return a.i == b.i && a.j == b.j; }

cluster_t *init_cluster(int radius) {
  int diameter = 2 * radius + 1;
  cluster_t *g = malloc(sizeof(cluster_t));
  if (g == NULL)
    return NULL;
  g->matrix = malloc(sizeof(bool*) * diameter);
  g->radius = radius;
  g->diameter = diameter;
  int i;
  for (i = 0; i < diameter; i++)
    g->matrix[i] = calloc(diameter, sizeof(bool));

  // start off with just the center point in the cluster
  g->matrix[radius][radius] = 1;

  return g;
}

// return coordinate of starting point for particle around circle of size radius
coord_t create_start(int rb, random_t *seed) {
  float deg = next_random_float(seed, TAU);
  coord_t res;
  res.i = (int)(rb * sin(deg));
  res.j = (int)(rb * cos(deg));
  return res;
}

coord_t create_start_large(int rb, random_t *seed) {
  float deg = omp_get_thread_num() * (TAU / omp_get_num_threads()) + next_random_float(seed, TAU / omp_get_num_threads());
  coord_t res;
  res.i = (int)(rb * sin(deg));
  res.j = (int)(rb * cos(deg));
  return res;
}

// create a walk of length k
void create_walk(coord_t *walk, int k, coord_t start, random_t *seed) {
  int i;
  for (i = 1; i < k; i++) {
    random_t f = next_random_float(seed, 4);
    int dir = round(f);
    switch (dir) {
    case 0:
      walk[i].i = 1;
      walk[i].j = 0;
      break;
    case 1:
      walk[i].i = 0;
      walk[i].j = 1;
      break;
    case 2:
      walk[i].i = -1;
      walk[i].j = 0;
      break;
    case 3:
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

// initialize the starting radius, number of walkers, and walk length
param_t *step_1(int rc, int M) {
  param_t *ret = malloc(sizeof(param_t));
  ret->rb = rc + 2 + M / 10;
  ret->k = choose_k(rc);
  ret->w = choose_w(M);
  return ret;
}

// create each walk
coord_t **step_2(param_t *params, random_t *seeds) {
  int i;
  int w = params->w;
  coord_t **all_walks = malloc(sizeof(coord_t *) * w);
#if OMP
#pragma omp parallel for
#endif
  for (i = 0; i < w; i++) {
    all_walks[i] = malloc(sizeof(coord_t) * params->k);
    create_walk(all_walks[i], params->k,
                create_start(params->rb, seeds + i), seeds + i);
  }
  return all_walks;
}

// create each walk
coord_t **step_2_large(param_t *params, random_t *seeds) {
  int i;
  int w = params->w;
  coord_t **all_walks = malloc(sizeof(coord_t *) * w);
  for (i = 0; i < w; i++) {
    all_walks[i] = malloc(sizeof(coord_t) * params->k);
    create_walk(all_walks[i], params->k,
                create_start_large(params->rb, seeds + i), seeds + i);
  }
  return all_walks;
}

// return true iff particle sticks to the cluster
bool is_sticky(coord_t loc, cluster_t *cluster) {
  int locx = loc.i + cluster->radius;
  int locy = loc.j + cluster->radius;
  bool **matrix = cluster->matrix;
  return
    locx >= 1 && locx < cluster->diameter - 1
    && locy >= 1 && locy < cluster->diameter - 1
    && (matrix[locx - 1][locy]
    ||  matrix[locx + 1][locy]
    ||  matrix[locx][locy - 1]
    ||  matrix[locx][locy + 1]
    ||  matrix[locx + 1][locy + 1]
    ||  matrix[locx + 1][locy - 1]
    ||  matrix[locx - 1][locy + 1]
    ||  matrix[locx - 1][locy - 1]
    );
}

// return the points where each particle sticks to the cluster (-1 if it doesn't stick)
int *step_3(coord_t **walks, cluster_t *cluster, param_t *params) {
  int i, j;
  int *res = malloc(sizeof(int) * params->w);
  memset(res, -1, params->w * sizeof(int));
#if OMP
#pragma omp parallel for
#endif
  for (i = 0; i < params->w; i++) {
    for (j = 0; j < params->k; j++) {
      if (is_sticky(walks[i][j], cluster)) {
        res[i] = j;
        break;
      }
    }
  }
  return res;
}

// return the first particle id where two particles interfere
// takes the higher of the two ids
// returns -1 if there is no interference
int step_4(int *res, coord_t **walks, param_t *params) {
  int i, j, k;
  for (i = 0; i < params->w; i++) {
    for (j = 0; j < i; j++) {
      if (res[j] == -1)
        break;
      for (k = 0; k < params->k; k++) {
        if (equal_coord(walks[j][res[j]], walks[i][k]))
          return i;
      }
    }
  }
  return -1;
}

// add particles that stick to the cluster
int step_5(cluster_t *cluster, int *res, coord_t **walks, param_t *params,
           int k, int *M) {
  int i;
  int rc = 0;
  // if k is -1 then every particle that sticks should be added to the cluster
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
    }
  }
  return cluster->radius < rc ? cluster->radius : rc;
}

void output_cluster(cluster_t *g, FILE *outfile) {
  int i, j;
  for (i = 0; i < g->diameter; i++) {
    for (j = 0; j < g->diameter; j++) {
      if (g->matrix[i][j] == 0)
        fprintf(outfile, "%d %d\n", i, j);
    }
  }
}

random_t *init_seeds(random_t seed, int max_radius) {
  int seeds_size = KFACTOR * KFACTOR * max_radius * max_radius;
  random_t *seeds = malloc(sizeof(random_t) * seeds_size);
  random_t *seedseed = malloc(sizeof(random_t));
  *seedseed = DEFAULTSEED;
  int i;
  for (i = 0; i < seeds_size; i++)
    seeds[i] = (random_t)next_random_float(seedseed, seeds_size);

  return seeds;
}

void do_simulation(cluster_t *cluster, random_t *seeds, int max_radius) {
  int rc = 1;
  int M = 1;

  while (rc < max_radius) {
    if (rc < RAD_THRESHOLD || true) {
      param_t *params = step_1(rc, M);
      START_ACTIVITY(ACTIVITY_CREATE_WALKS);
      coord_t **walks = step_2(params, seeds);
      FINISH_ACTIVITY(ACTIVITY_CREATE_WALKS);

      START_ACTIVITY(ACTIVITY_FIND_STICKING);
      int *res = step_3(walks, cluster, params);
      FINISH_ACTIVITY(ACTIVITY_FIND_STICKING);

      START_ACTIVITY(ACTIVITY_FIND_INTERFERENCE);
      int k = step_4(res, walks, params);
      FINISH_ACTIVITY(ACTIVITY_FIND_INTERFERENCE);

      START_ACTIVITY(ACTIVITY_ADD_TO_CLUSTER);
      int thisrc = step_5(cluster, res, walks, params, k, &M);
      FINISH_ACTIVITY(ACTIVITY_ADD_TO_CLUSTER);

      rc = thisrc > rc ? thisrc : rc;
      int i;
      for (i = 0; i < params->w; i++)
        free(walks[i]);
      free(walks);
      free(res);
      free(params);
    } else {
#if OMP
#pragma omp parallel
#endif
      {
        random_t *seeds = init_seeds(DEFAULTSEED, max_radius);
        param_t *params = step_1(rc, M);
        params->w /= omp_get_num_threads();
        coord_t **walks = step_2_large(params, seeds);
        int *res = step_3(walks, cluster, params);
        int k = step_4(res, walks, params);
        int thisrc = step_5(cluster, res, walks, params, k, &M);
        rc = thisrc > rc ? thisrc : rc;
        int i;
        for (i = 0; i < params->w; i++)
          free(walks[i]);
        free(walks);
        free(res);
        free(params);
        free(seeds);
      }
    }
  }
}

void simulate(bool tracking, char* image_name, char *test_output_name, int max_radius, random_t seed) {
  track_activity(tracking);
  cluster_t *cluster = init_cluster(max_radius);
  random_t *seeds = init_seeds(seed, max_radius);

  do_simulation(cluster, seeds, max_radius);

  if (image_name != NULL)
    generate_image(cluster->matrix, cluster->diameter, cluster->radius, image_name);

  if (test_output_name != NULL) {
    FILE *outfile = fopen(test_output_name, "a+");
    output_cluster(cluster, outfile);
    fclose(outfile);
  }


  SHOW_ACTIVITY(stderr, tracking);
}
