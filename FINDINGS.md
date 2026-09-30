# Findings

**Setup:** Intel Core i9-13900K (32 threads), NVIDIA GeForce RTX 4090 (24 GB), Linux 6.8, OpenJDK 21.0.2,
TornadoVM 7.0.1 (CUDA backend), BioJava 7.3.0. Measured 2026-09-30. Times are wall clock around the public
call (host preparation, host-to-device and device-to-host transfers, and kernel), taking the best of 3 runs
after one warm-up. The one-off kernel JIT (about 0.5 s per JVM) is excluded and reported separately.

**Two CPU baselines**, so the GPU number is not inflated by BioJava's implementation:
* **BioJava**: the stock code, multi-threaded (`AsaCalculator` with 32 threads; `Alignments` runs the
  pairs on its 32-thread pool).
* **Lean CPU**: the *same algorithm as the GPU kernel* in plain Java on primitive arrays, one task per
  atom or pair on all 32 cores (`asa/CpuAsa`, `align/CpuAlignmentScores`; the library's own CPU path). This is the fair GPU-vs-CPU
  comparison. BioJava vs lean CPU shows how much BioJava leaves on the table in pure Java.

**Correctness:** every run below checks the GPU result against BioJava's (the alignment rows also include
linear gap penalties in the tests): per-atom ASA is bit-identical
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
| 40 seqs, 100-300 aa | 780 | 0.03 G | GLOBAL | 242 ms | 12.4 ms | 5.2 ms | 46x | 2.4x | 6 |
| 40 seqs, 100-300 aa | 780 | 0.03 G | LOCAL | 245 ms | 6.8 ms | 6.2 ms | 40x | 1.1x | 5 |
| 200 seqs, 100-500 aa | 19,900 | 1.84 G | GLOBAL | 9.72 s | 381 ms | 17.2 ms | 567x | 22x | 107 |
| 200 seqs, 100-500 aa | 19,900 | 1.84 G | LOCAL | 10.2 s | 257 ms | 16.7 ms | 615x | 15x | 110 |
| Pfam PF00104 | 39,903 | 1.02 G | GLOBAL | 6.04 s | 246 ms | 8.9 ms | 680x | 28x | 115 |
| Pfam PF00104 | 39,903 | 1.02 G | LOCAL | 7.02 s | 199 ms | 10.0 ms | 705x | 20x | 102 |
| 1000 seqs, 200-500 aa | 499,500 | 61.8 G | GLOBAL | not run | 12.7 s | 260 ms | - | 49x | 238 |
| 1000 seqs, 200-500 aa | 499,500 | 61.8 G | LOCAL | not run | 9.1 s | 261 ms | - | 35x | 237 |
| 2000 seqs, 100-400 aa | 1,999,000 | 126 G | GLOBAL | not run | 34.1 s | 608 ms | - | 56x | 208 |
| 2000 seqs, 100-400 aa | 1,999,000 | 126 G | LOCAL | not run | 24.0 s | 609 ms | - | 39x | 207 |

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
* The execution plan is persistent and capacity-padded, as for ASA, so a call costs ~4 ms of fixed overhead.
  The GPU wins from about 1k pairs of 200 aa (0.03 G cells); below the 1e7-cell threshold the library runs
  the lean CPU code (e.g. 45 pairs: CPU 1.5 ms vs GPU 3.9 ms).
* History: the first version rebuilt its plan per call (~20 ms) and had an overflowing pair-sort key, which
  scrambled the warp ordering; it measured 4-34x vs lean CPU. Fixing both gave the table above.
* Guide trees built from the GPU scorers are identical (Newick string) to BioJava's (test on 60 Pfam
  sequences).

## External baselines (C libraries, same inputs)

`scripts/external/compare.sh` builds the libraries and runs them on exactly the inputs of our benchmarks.

**parasail** (SIMD striped NW/SW with 8/16/32-bit saturation-checked profiles, OpenMP, 32 threads, `-O3
-march=native`). BioJava's gaps (open 10, extend 1: a gap of k costs 10 + k) are parasail's open 11, extend 1.
**All scores are identical** across BioJava, our GPU kernel and parasail.

| Set | Type | parasail 32 thr | Lean Java 32 thr | GPU (CUDA) | GPU vs parasail |
|---|---|---:|---:|---:|---:|
| 200 seqs, 19,900 pairs | NW | 52.0 ms | 381 ms | 17.1 ms | 3.0x |
| 200 seqs, 19,900 pairs | SW | 11.5 ms | 257 ms | 16.9 ms | 0.7x |
| Pfam PF00104, 39,903 pairs | NW | 61.9 ms | 246 ms | 8.6 ms | 7.2x |
| Pfam PF00104, 39,903 pairs | SW | 18.8 ms | 199 ms | 8.9 ms | 2.1x |
| 1000 seqs, 499,500 pairs | NW | 1580 ms | 12.7 s | 259 ms | 6.1x |
| 1000 seqs, 499,500 pairs | SW | 261 ms | 9.1 s | 259 ms | 1.0x |
| 2000 seqs, 1,999,000 pairs | NW | 3996 ms | 34.1 s | 598 ms | 6.7x |
| 2000 seqs, 1,999,000 pairs | SW | 647 ms | 24.0 s | 599 ms | 1.1x |

* Against a state-of-the-art SIMD CPU library the GPU kernel is **3-7x faster for global alignment** and only
  **0.7-2x for local alignment**. parasail's local alignment of unrelated sequences stays in 8-bit lanes (32
  per AVX2 vector), which our 32-bit-per-cell kernel does not match.
* Our kernel is at ~230-280 GCUPS. Packed 8/16-bit arithmetic on the GPU (as CUDASW++ does) is the known
  route to more; not attempted. Kernel variants tried: strip heights 8/16/32 (32 best); a fully unrolled strip
  with the row state in registers instead of local memory was 10-20% slower; batches sized to 4 GB of
  row buffers (more resident threads) gave +10%. The profiler shows the kernel is ~85% of the time.

**FreeSASA** (C, Shrake-Rupley, 1000 points, probe 1.4, called through its C API with BioJava's atoms and
radii; its S&R supports at most 16 threads). Totals agree with BioJava's within 0.05% (FreeSASA uses its own
sphere points).

| PDB | Atoms | FreeSASA 16 thr | Lean Java 16 thr | Lean Java 32 thr | GPU | GPU vs FreeSASA |
|---|---:|---:|---:|---:|---:|---:|
| 1SMT | 1,574 | 7.0 ms | - | 1.7 ms | 0.8 ms | 8.8x |
| 4HHB | 4,384 | 13.3 ms | - | 2.8 ms | 1.7 ms | 7.8x |
| 1AON | 58,674 | 144 ms | 39.7 ms | 28.7 ms | 12.3 ms | 11.7x |
| 4V6X | 237,685 | 561 ms | 175 ms | 109 ms | 39.2 ms | 14.3x |
| 3J3Q | 2,440,800 | 5.29 s | 1.64 s | 1.13 s | 372 ms | 14.2x |

* At equal thread counts the lean Java ASA is ~3.2x faster than FreeSASA (neighbours sorted by distance, so
  the occlusion test usually stops at the first neighbour).

## Backends and devices

The same code and the same 20 parity tests, on the other TornadoVM backends available on this machine (OpenCL SDK
`tornadovm-7.0.1-jdk21-opencl`; device picked with `-Dsw.scores.device=0:1` / `-Dasa.fused.device=...`):

| Device (backend) | Parity tests | Alignment, 20k pairs (1.84 G cells) | Alignment, 500k pairs (61.8 G cells) | ASA |
|---|---|---|---|---|
| RTX 4090 (CUDA/PTX) | 20/20 exact | 17.2 ms, 107 GCUPS | 260 ms, 238 GCUPS | as above |
| RTX 4090 (OpenCL) | 20/20 exact | 15.0 ms, 122 GCUPS | 230 ms, 269 GCUPS (65x lean CPU) | runs, exact |
| Intel UHD 770 iGPU (OpenCL) | 20/20 exact | 3.28 s, 0.6 GCUPS | not run | no FP64: falls back to the CPU path |

* On the 4090, OpenCL is 10-15% faster than CUDA/PTX for this integer kernel.
* The Intel iGPU gives exact scores but is 5-10x slower than the 32-core CPU: not worth using. It lacks FP64,
  so the ASA kernel is refused at compile time (`TornadoDeviceFP64NotSupported`); the library then runs ASA on
  the CPU path, and keeps the alignment kernel on the iGPU (the fallback is per operation).

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
* Identity/similarity scorer types (`*_IDENTITIES`, the MSA default) need a traceback on the GPU.
* AMD and Apple GPUs are not measured yet. FP64-less devices would need an FP32 ASA kernel with an exact
  re-check of the boundary cases.
* JMH harness: the benchmarks are simple best-of-N timers.
