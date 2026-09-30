# Findings (work in progress: final tables being regenerated)

Machine: Intel Core i9-13900K (32 threads), NVIDIA RTX 4090, JDK 21.0.2, TornadoVM 7.0.1 (CUDA), BioJava 7.3.0.

Every speedup is reported against **two** CPU baselines:
* **BioJava**: the stock BioJava code, multi-threaded (`AsaCalculator` with 32 threads; `Alignments` thread pool).
* **lean CPU**: the *same algorithm as the GPU kernel*, in plain Java on primitive arrays, on all 32 cores
  (`bench/AsaLeanCpu`, `bench/CpuScores`). This is the fair GPU-vs-CPU comparison; BioJava vs lean CPU
  measures how much BioJava leaves on the table in pure Java.

All results are exact: bit-identical per-atom ASA, identical integer alignment scores.

## CE (jCE): dropped
`CECalculator.initSumOfDistances` was ported and gives identical alignments, but it is only ~0-2% of CE's
runtime, so the end-to-end speedup was 0.98-1.01x. A JFR profile of CeMain (1CDG.A vs 1CGT.A) puts ~71% of
the samples in `CECalculator.dpAlign` and ~26% in `getScoreFromDistanceMatrices`: those are the real targets.

## CUDA libraries (cuBLAS, CUTLASS, cuFFT, cuSPARSE, cuDNN, cuRAND, cuDF)
None fits these kernels (geometric neighbour search with early exit; integer DP recurrences). A cuBLAS fit in
BioJava would be all-vs-all superposition (QCP cross-covariances of many models as one GEMM); not attempted.
