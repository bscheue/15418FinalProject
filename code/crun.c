#include <getopt.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stdio.h>

#include "rutil.h"


static void usage(char *name) {
  char *use_string = "[-m MAX_RADIUS] [-o OUTPUT_NAME] [-r IMAGE_NAME] [-s SEED] [-t]";
  printf("Usage: %s %s\n", name, use_string);
}

void simulate(bool tracking, char* image_name, char *test_output_name, int max_radius, random_t seed);

int main(int argc, char *argv[]) {
  int c;
  char *optstring = "hm:r:s:t";
  random_t seed = DEFAULTSEED;
  int max_radius = 10;
  bool tracking = false;
  char* image_name = NULL;
  char* test_output_name = NULL;
  while ((c = getopt(argc, argv, optstring)) != -1) {
    switch (c) {
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
      tracking = true;
      break;
    default:
      printf("Unknown option '%c'\n", c);
      usage(argv[0]);
      exit(1);
    }
  }

  simulate(tracking, image_name, test_output_name, max_radius, seed);
  return 0;
}
