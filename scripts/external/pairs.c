// All-pairs NW/SW scores with parasail (striped SIMD, saturation-checked profiles), OpenMP over queries.
// usage: pairs <fasta> <nw|sw> <open> <extend> <out-scores>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <omp.h>
#include "parasail.h"
#include "parasail/matrices/blosum62.h"

int main(int argc, char **argv) {
  FILE *f = fopen(argv[1], "r");
  int local = strcmp(argv[2], "sw") == 0, open = atoi(argv[3]), ext = atoi(argv[4]);
  int cap = 1024, n = 0; char **seq = malloc(cap * sizeof(char *)); int *len = malloc(cap * sizeof(int));
  char line[1 << 16]; char *cur = NULL; size_t curLen = 0;
  while (fgets(line, sizeof line, f)) {
    line[strcspn(line, "\r\n")] = 0;
    if (line[0] == '>') {
      if (n == cap) { cap *= 2; seq = realloc(seq, cap * sizeof(char *)); len = realloc(len, cap * sizeof(int)); }
      seq[n] = calloc(1, 1); len[n] = 0; n++;
    } else {
      size_t l = strlen(line);
      seq[n - 1] = realloc(seq[n - 1], len[n - 1] + l + 1);
      memcpy(seq[n - 1] + len[n - 1], line, l + 1); len[n - 1] += l;
    }
  }
  long nPairs = (long) n * (n - 1) / 2;
  int *scores = malloc(nPairs * sizeof(int));
  long *start = malloc(n * sizeof(long)); start[0] = 0;
  for (int i = 1; i < n; i++) start[i] = start[i - 1] + (n - i);
  double t0 = omp_get_wtime();
#pragma omp parallel for schedule(dynamic, 1)
  for (int i = 0; i < n - 1; i++) {
    parasail_profile_t *prof = parasail_profile_create_sat(seq[i], len[i], &parasail_blosum62);
    for (int j = i + 1; j < n; j++) {
      parasail_result_t *r = local ? parasail_sw_striped_profile_sat(prof, seq[j], len[j], open, ext)
                                   : parasail_nw_striped_profile_sat(prof, seq[j], len[j], open, ext);
      scores[start[i] + j - i - 1] = parasail_result_get_score(r);
      parasail_result_free(r);
    }
    parasail_profile_free(prof);
  }
  double t1 = omp_get_wtime();
  long cells = 0;
  for (int i = 0; i < n; i++) for (int j = i + 1; j < n; j++) cells += (long) len[i] * len[j];
  printf("parasail %s: %ld pairs, %.2f Gcells, %.1f ms, %.1f GCUPS, %d threads\n", argv[2], nPairs, cells / 1e9,
         (t1 - t0) * 1e3, cells / (t1 - t0) / 1e9, omp_get_max_threads());
  FILE *o = fopen(argv[5], "w");
  for (long p = 0; p < nPairs; p++) fprintf(o, "%d\n", scores[p]);
  fclose(o);
  return 0;
}
