# Findings

**Setup:** Intel Core i9-13900K (32 threads), NVIDIA GeForce RTX 4090 (24 GB), Linux 6.8, OpenJDK 21.0.2,
TornadoVM 7.0.1 (CUDA backend), BioJava 7.3.0. Measured 2026-09-30. Times are wall clock around the public
call (host preparation, host-to-device and device-to-host transfers, and kernel), taking the best of 3 runs
after one warm-up. The one-off kernel JIT (about 0.5 s per JVM) is excluded and reported separately.

**Two CPU baselines**, so the GPU number is not inflated by BioJava's implementation:
* **BioJava**: the stock code, multi-threaded (`AsaCalculator` with 32 threads; `Alignments` runs the
  pairs on its 32-thread pool).
* **Lean CPU**: the *same algorithm as the GPU kernel* in plain Java on primitive arrays, one task per
  atom or pair on all 32 cores (`bench/AsaLeanCpu`, `bench/CpuScores`). This is the fair GPU-vs-CPU
  comparison. BioJava vs lean CPU shows how much BioJava leaves on the table in pure Java.

**Correctness:** every run below checks the GPU result against BioJava's: per-atom ASA is bit-identical
(0 differing atoms) and alignment scores are identical integers. The lean CPU baselines are also exact.
The 16 JUnit parity tests pass under TornadoVM (`scripts/test-gpu.sh`) and on a plain JVM (fallback).

## ASA (Shrake-Rupley, 1000 sphere points, probe 1.4 A)

`bench/AsaBench`, all non-H atoms of model 1.

| PDB | Atoms | BioJava 1 thread | BioJava 32 threads | Lean CPU 32 | GPU | vs BioJava 32 | vs lean CPU |
|---|---:|---:|---:|---:|---:|---:|---:|
| 1SMT | 1,574 | 13.1 ms | 6.5 ms | 1.7 ms | 0.8 ms | 7.9x | 2.1x |
| 4HHB | 4,384 | 30.1 ms | 12.7 ms | 2.8 ms | 1.7 ms | 7.4x | 1.6x |
| 1AON | 58,674 | 364 ms | 132 ms | 27.6 ms | 12.2 ms | 10.8x | 2.3x |
| 4V6X | 237,685 | 1.47 s | 519 ms | 119 ms | 35.3 ms | 14.7x | 3.4x |
| 3J3Q | 2,440,800 | 16.6 s | 6.8 s | 1.19 s | 390 ms | 17.5x | 3.1x |

* Kernel design: one work-group (128 threads) per atom; cell-grid neighbour search, compacted with a
  work-group prefix sum; rank sort by distance; early-exit occlusion test on the local-memory neighbour
  list. Everything is in double precision with BioJava's operation order, hence bit-exact.
* The execution plan is persistent and capacity-padded: compiled once, re-executed for every structure.
  Before that change each call rebuilt its plan and paid about 30 ms. Re-executing a plan costs 0.6 ms.
* TornadoVM profiler on 3J3Q: kernel 272 ms, copy-in 8.5 ms; the remainder is host-side cell bucketing.
* Tried and dropped: an FP32 occlusion test with FP64 re-check near the boundary. It stayed exact but gave
  no speedup, because the kernel is not FP64-throughput bound.
* Side finding: the lean pure-Java version is **4-6x faster than BioJava's multi-threaded `AsaCalculator`**
  with identical results. That is a CPU-only improvement BioJava could adopt independently of GPUs.

## All-pairs alignment scores (Needleman-Wunsch `GLOBAL`, Smith-Waterman `LOCAL`)

`bench/SwBench`, BLOSUM62, gap open 10 / extend 1 (BioJava defaults). Random protein sets use a fixed seed and
a uniform length range; Pfam PF00104 is 283 ungapped family members (about 160 aa on average: 1.02 G cells over 39,903 pairs).

| Set | Pairs | DP cells | Type | BioJava (32 thr) | Lean CPU 32 | GPU | vs BioJava | vs lean CPU | GPU GCUPS |
|---|---:|---:|---|---:|---:|---:|---:|---:|---:|
| 10 seqs, 100-300 aa | 45 | < 0.01 G | GLOBAL | 44 ms | 2.0 ms | 28 ms | 1.6x | 0.07x | 0.1 |
| 40 seqs, 100-300 aa | 780 | 0.03 G | GLOBAL | 207 ms | 10.8 ms | 22 ms | 9.3x | 0.5x | 1.5 |
| 200 seqs, 100-500 aa | 19,900 | 1.84 G | GLOBAL | 9.65 s | 379 ms | 59 ms | 163x | 6.4x | 31 |
| 200 seqs, 100-500 aa | 19,900 | 1.84 G | LOCAL | 10.1 s | 285 ms | 58 ms | 174x | 4.9x | 32 |
| Pfam PF00104 | 39,903 | 1.02 G | GLOBAL | 5.60 s | 205 ms | 48 ms | 117x | 4.3x | 21 |
| Pfam PF00104 | 39,903 | 1.02 G | LOCAL | 6.59 s | 228 ms | 48 ms | 138x | 4.8x | 21 |
| 1000 seqs, 200-500 aa | 499,500 | 61.8 G | GLOBAL | not run | 14.6 s | 616 ms | - | 23.6x | 100 |
| 1000 seqs, 200-500 aa | 499,500 | 61.8 G | LOCAL | not run | 11.2 s | 557 ms | - | 20.1x | 111 |
| 2000 seqs, 100-400 aa | 1,999,000 | 126 G | GLOBAL | not run | 30.0 s | 877 ms | - | 34.2x | 144 |
| 2000 seqs, 100-400 aa | 1,999,000 | 126 G | LOCAL | not run | 22.2 s | 896 ms | - | 24.8x | 141 |

* BioJava was not run on the largest sets. At its measured ~0.19 GCUPS it would need an estimated
  5-11 minutes each; that is an extrapolation, not a measurement.
* BioJava's aligner keeps a full traceback (a `Last[]` object per DP cell) even when only the score is
  wanted, which is why it is about 25x slower than the lean CPU code. Most of the "vs BioJava" factor is
  therefore algorithmic (score-only, primitive arrays), not the GPU.
* The lean CPU baseline varies by about ±20% between runs (turbo/thermal); the table shows one run's best of 3.
* Kernel design: one pair per thread; pairs sorted by DP size so a warp does similar work; the DP boundary
  row lives in interleaved (coalesced) global buffers; the DP is strip-mined over 32 query rows with the
  per-row state in local memory. A naive one-row-at-a-time version reached only 7.5-9.6 GCUPS
  (0.8-2.2x the lean CPU). Strip height 8 gave 62-72 GCUPS, 16 gave 78-98 and 32 gave 83-144.
* Small batches lose to the CPU: the fixed cost is about 20 ms per call, because the plan is currently rebuilt
  per call and the batch is below the size that fills the GPU. Hence the 5e7-cell threshold.
* Guide trees built from the GPU scorers are identical (Newick string) to BioJava's (test on 60 Pfam
  sequences).

## Investigated and not included

**CE structure alignment (jCE).** `CECalculator.initSumOfDistances` was ported; it gave the same
alignments and scores (the distance matrix differed in the last bits, likely FMA contraction), but it is only about 0-2% of CE's runtime, so the end-to-end speedup was 0.98-1.01x (1CDG/1CGT,
1TIM/1CDG, 5A22/6U1X). A JFR profile of `CeMain` (1CDG.A vs 1CGT.A) puts about 71% of samples in
`CECalculator.dpAlign` and about 26% in `getScoreFromDistanceMatrices`: those are the real targets.

**TornadoVM CUDA library bindings** (cuBLAS/cuBLASLt, CUTLASS, cuFFT, cuSPARSE, cuDNN, cuRAND, cuDF). None
fits these kernels: ASA is geometric neighbour search with early exit, alignment is an integer DP recurrence,
and cuRAND would change BioJava's deterministic spiral points. A natural cuBLAS use in BioJava would be
all-vs-all superposition: the QCP cross-covariance matrices of many models form one (3M x N)(N x 3M) GEMM.
Not attempted.

## Not yet done
* Persistent (cached) execution plan for the alignment path, as for ASA, to cut the ~20 ms per-call cost.
* Identity/similarity scorer types (`*_IDENTITIES`, the MSA default) need a traceback on the GPU.
* Other GPUs and backends (OpenCL, AMD, Intel, Apple) are not measured yet.
* JMH harness: the benchmarks are simple best-of-N timers.
