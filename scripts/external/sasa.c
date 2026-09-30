// FreeSASA Shrake-Rupley on BioJava's atoms and radii. usage: sasa <xyzr-file> <threads> <reps>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "freesasa.h"
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec * 1e-9; }
int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r"); int n; if (fscanf(f, "%d", &n) != 1) return 1;
  double *xyz = malloc(3 * n * sizeof(double)), *r = malloc(n * sizeof(double));
  for (int i = 0; i < n; i++) if (fscanf(f, "%lf %lf %lf %lf", &xyz[3*i], &xyz[3*i+1], &xyz[3*i+2], &r[i]) != 4) return 2;
  freesasa_parameters p = freesasa_default_parameters;
  p.alg = FREESASA_SHRAKE_RUPLEY; p.shrake_rupley_n_points = 1000; p.probe_radius = 1.4; p.n_threads = atoi(argv[2]);
  double best = 1e30, total = 0;
  for (int k = 0; k < atoi(argv[3]); k++) {
    double t0 = now();
    freesasa_result *res = freesasa_calc_coord(xyz, r, n, &p);
    double t = now() - t0; if (t < best) best = t;
    total = res->total; freesasa_result_free(res);
  }
  printf("freesasa SR1000 %d thr: %d atoms, %.1f ms, total %.2f\n", p.n_threads, n, best * 1e3, total);
  return 0;
}
