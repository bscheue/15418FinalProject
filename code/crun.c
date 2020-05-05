#include <getopt.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

#include "rutil.h"

#if OMP
#include <omp.h>
#else
#include "fake_omp.h"
#endif


static void usage(char *name) {
  char *use_string = "[-b] [-m MAX_RADIUS] [-o TEST_OUTPUT_NAME] [-r IMAGE_NAME] [-s SEED] [-t NUM_THREADS]";
  printf("Usage: %s %s\n", name, use_string);
}

void simulate(bool tracking, char* image_name, char *test_output_name, int max_radius, random_t seed);

int main(int argc, char *argv[]) {
  int c;
  char *optstring = "bhm:o:r:s:t:";
  random_t seed = DEFAULTSEED;
  int max_radius = 125;
  bool tracking = false;
  char* image_name = NULL;
  char* test_output_name = NULL;
  int thread_count = 1;
  while ((c = getopt(argc, argv, optstring)) != -1) {
    switch (c) {
    case 'b':
      tracking = true;
      break;
    case 'h':
      usage(argv[0]);
      break;
    case 'm':
      max_radius = atoi(optarg);
      break;
    case 'o':
      test_output_name = optarg;
      break;
    case 'r':
      image_name = optarg;
      break;
    case 's':
      seed = atoi(optarg);
      break;
    case 't':
      thread_count = atoi(optarg);
      break;
    default:
      printf("Unknown option '%c'\n", c);
      usage(argv[0]);
      exit(1);
    }
  }

  omp_set_num_threads(thread_count);

  simulate(tracking, image_name, test_output_name, max_radius, seed);
  return 0;
}
