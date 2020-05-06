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

#define EMPTY 0
#define STICKY 1
#define OCCUPIED 2

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

    char **matrix; // matrix

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
  g->matrix = malloc(sizeof(char*) * diameter);
  g->radius = radius;
  g->diameter = diameter;
  int i;
  for (i = 0; i < diameter; i++)
    g->matrix[i] = calloc(diameter, sizeof(char));

  // start off with just the center point in the cluster
  g->matrix[radius][radius] = OCCUPIED;

  g->matrix[radius-1][radius] = STICKY;
  g->matrix[radius-1][radius-1] = STICKY;
  g->matrix[radius-1][radius+1] = STICKY;
  g->matrix[radius+1][radius] = STICKY;
  g->matrix[radius+1][radius-1] = STICKY;
  g->matrix[radius+1][radius+1] = STICKY;
  g->matrix[radius][radius+1] = STICKY;
  g->matrix[radius][radius-1] = STICKY;

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
// create a walk of length k
void create_walk(coord_t *walk, int k, coord_t start, random_t *seed) {
  int i;
  walk[0].i = start.i;
  walk[0].j = start.j;
  int accum_i = start.i;
  int accum_j = start.j;
  for (i = 1; i < k; i++) {
    random_t f = next_random_float(seed, 4);
    int dir = round(f);
    switch (dir) {
    case 0:
      walk[i].i = accum_i + 1;
      walk[i].j = accum_j;
      break;
    case 1:
      walk[i].j = accum_j + 1;
      walk[i].i = accum_i;
      break;
    case 2:
      walk[i].i = accum_i - 1;
      walk[i].j = accum_j;
      break;
    case 3:
      walk[i].j = accum_j - 1;
      walk[i].i = accum_i;
      break;
    default:
      printf("mistake\n");
    }
    accum_i = walk[i].i;
    accum_j = walk[i].j;
  }
}

// initialize the starting radius, number of walkers, and walk length
param_t *step_1(param_t *params, int rc, int M) {
  params->rb = rc + 2 + M / 10;
  params->k = choose_k(rc);
  params->w = choose_w(M);
  return params;
}

// create each walk
void step_2(coord_t **walks, param_t *params, random_t *seeds) {
  int i;
#if OMP
#pragma omp parallel for
#endif
  for (i = 0; i < params->w; i++) {
    create_walk(walks[i], params->k,
                create_start(params->rb, seeds + i), seeds + i);
  }
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
  char **matrix = cluster->matrix;
  if (locx <= 0 || locx >= cluster->diameter-1
    || locy <= 0 || locy >= cluster->diameter-1) return false;
  return matrix[locx][locy] == STICKY;
}

// return the points where each particle sticks to the cluster (-1 if it doesn't stick)
void step_3(int *res, coord_t **walks, cluster_t *cluster, param_t *params) {
  int i, j;
  memset(res, -1, params->w * sizeof(int));
#if OMP
#pragma omp parallel for schedule(guided)
#endif
  for (i = 0; i < params->w; i++) {
    for (j = 0; j < params->k; j++) {
      if (is_sticky(walks[i][j], cluster)) {
        res[i] = j;
        break;
      }
    }
  }
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

bool add_particle(cluster_t *cluster, int i, int j) {
  int locx = i + cluster->radius;
  int locy = j + cluster->radius;
  char **matrix = cluster->matrix;

  if (matrix[locx][locy] == OCCUPIED
    || locx <= 0 || locx >= cluster->diameter - 1
    || locy <= 0 || locy >= cluster->diameter - 1) {
    return false;
  }

  matrix[locx][locy] = OCCUPIED;
  if (locx >= 1) {
    if (matrix[locx-1][locy] == EMPTY) matrix[locx-1][locy] = STICKY;
    if (locy >= 1 && matrix[locx-1][locy-1] == EMPTY) matrix[locx-1][locy-1] = STICKY;
    if (locy < cluster->diameter - 1 && matrix[locx-1][locy+1] == EMPTY) matrix[locx-1][locy+1] = STICKY;
  }
  if (locx < cluster->diameter-1) {
    if (matrix[locx+1][locy] == EMPTY) matrix[locx+1][locy] = STICKY;
    if (locy >= 1 && matrix[locx+1][locy-1] == EMPTY) matrix[locx+1][locy-1] = STICKY;
    if (locy < cluster->diameter - 1 && matrix[locx+1][locy+1] == EMPTY) matrix[locx+1][locy+1] = STICKY;
  }
  if (locy >= 1 && matrix[locx][locy-1] == EMPTY) matrix[locx][locy-1] = STICKY;
  if (locy < cluster->diameter-1 && matrix[locx][locy+1] == EMPTY) matrix[locx][locy+1] = STICKY;
  return true;

}

// add particles that stick to the cluster
int step_5(cluster_t *cluster, int *res, coord_t **walks, param_t *params,
           int k, int *M) {
  int i;
  int rc = 0;
  // if k is -1 then every particle that sticks should be added to the cluster
  k = k == -1 ? params->w : k;
#if OMP
#pragma omp parallel for
#endif 
  for (i = 0; i < k; i++) {
    if (res[i] == -1) {
      continue;
    }
    coord_t tup = walks[i][res[i]];
    int this_rc = round(sqrt(tup.i * tup.i + tup.j * tup.j));
#if OMP
#pragma omp critical 
#endif 
{
    if (add_particle(cluster, tup.i, tup.j)) {
      rc = this_rc > rc ? this_rc : rc;
      *M += 1;
    }
  }
}
  return cluster->radius < rc ? cluster->radius : rc;
}

// add particles that stick to the cluster
int step_5_large(cluster_t *cluster, int *res, coord_t **walks, param_t *params,
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
    if (add_particle(cluster, tup.i, tup.j) == true) {
#if OMP
#pragma omp critical
#endif
      rc = this_rc > rc ? this_rc : rc;
#if OMP
#pragma omp critical
#endif
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

  param_t *params = malloc(sizeof(param_t));
  int max_w = choose_w(cluster->diameter * cluster->diameter);
  int *res = malloc(sizeof(int) * max_w);
  coord_t **walks = malloc(sizeof(coord_t *) * max_w);
  int i;
  for (i = 0; i < max_w; i++)
    walks[i] = malloc(sizeof(coord_t) * choose_k(max_radius));

  while (rc < max_radius) {
    if (rc < RAD_THRESHOLD || true) {
      params = step_1(params, rc, M);
      START_ACTIVITY(ACTIVITY_CREATE_WALKS);
      step_2(walks, params, seeds);
      FINISH_ACTIVITY(ACTIVITY_CREATE_WALKS);

      START_ACTIVITY(ACTIVITY_FIND_STICKING);
      step_3(res, walks, cluster, params);
      FINISH_ACTIVITY(ACTIVITY_FIND_STICKING);

      START_ACTIVITY(ACTIVITY_FIND_INTERFERENCE);
      int k = step_4(res, walks, params);
      FINISH_ACTIVITY(ACTIVITY_FIND_INTERFERENCE);

      START_ACTIVITY(ACTIVITY_ADD_TO_CLUSTER);
      int thisrc = step_5(cluster, res, walks, params, k, &M);
      FINISH_ACTIVITY(ACTIVITY_ADD_TO_CLUSTER);
      rc = thisrc > rc ? thisrc : rc;
    }
    else
#if OMP
#pragma omp parallel
#endif
      {
        random_t *seeds = init_seeds(DEFAULTSEED, max_radius);
        param_t *params = step_1(params, rc, M);
        params->w /= omp_get_num_threads();
        coord_t **walks = step_2_large(params, seeds);
        step_3(res, walks, cluster, params);
        int k = step_4(res, walks, params);
        int thisrc = step_5_large(cluster, res, walks, params, k, &M);
        rc = thisrc > rc ? thisrc : rc;
        int i;
        for (i = 0; i < params->w; i++)
          free(walks[i]);
        free(walks);
        free(seeds);
      }
    }

  for (i = 0; i < max_w; i++)
    free(walks[i]);
  free(walks);
  free(res);
  free(params);
  free(seeds);
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
