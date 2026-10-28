# biojava-tornadovm-extensions

GPU acceleration of two [BioJava](https://github.com/biojava/biojava) hot spots with
[TornadoVM](https://github.com/beehive-lab/TornadoVM), as **drop-in classes that return exactly the same
results** as the BioJava code they replace:

| BioJava class | Drop-in replacement | Result |
|---|---|---|
| `org.biojava.nbio.structure.asa.AsaCalculator` | `org.biojava.tornado.asa.TornadoAsaCalculator` | bit-identical per-atom and per-residue ASA |
| `org.biojava.nbio.alignment.Alignments.getAllPairsScores` / `getAllPairsScorers` (`GLOBAL`, `LOCAL`) | `org.biojava.tornado.align.TornadoAlignments` | identical scores, identical `GuideTree` |

## Results at a glance

RTX 4090 vs a 32-thread i9-13900K, end-to-end (host↔device transfers included). All results are identical to
BioJava's. "Lean Java" is the same algorithm as the GPU kernel in plain Java on all cores; it is also this
library's CPU fallback. Full tables and methodology are in [FINDINGS.md](FINDINGS.md).

**All-pairs alignment scores** (BLOSUM62, gaps 10/1)

| Workload | BioJava (32 thr) | Lean Java (32 thr) | [parasail](https://github.com/jeffdaily/parasail), C SIMD (32 thr) | GPU |
|---|---:|---:|---:|---:|
| Global, 200 proteins, 20k pairs | 9.72 s | 381 ms | 52 ms | **17 ms** |
| Global, Pfam PF00104, 40k pairs | 6.04 s | 246 ms | 62 ms | **9 ms** |
| Global, 2,000 proteins, 2M pairs | ~11 min (est.) | 34.1 s | 4.0 s | **0.60 s** |
| Local, 2,000 proteins, 2M pairs | ~11 min (est.) | 24.0 s | 0.65 s | **0.60 s** |

**Accessible surface area** (Shrake-Rupley, 1000 points)

| Workload | BioJava (32 thr) | Lean Java (32 thr) | [FreeSASA](https://freesasa.github.io/), C (16 thr, its max) | GPU |
|---|---:|---:|---:|---:|
| 4V6X, 238k atoms | 519 ms | 109 ms | 561 ms | **39 ms** |
| 3J3Q, 2.4M atoms | 6.8 s | 1.13 s | 5.29 s | **0.37 s** |

The short version:
* The GPU is far ahead of BioJava today, and 3-7x ahead of parasail for global alignment.
* For local alignment the GPU only matches parasail.
* Without any GPU, the lean Java path alone is 4-6x (ASA) and ~25x (alignment) faster than BioJava.

**Tested devices:** RTX 4090 with the CUDA and OpenCL backends, and an Intel UHD 770 iGPU with OpenCL (no
FP64, so ASA uses the exact single-precision path). All tests give exact results on all three.

BioJava is not modified. This library is built against BioJava 7.3.0 from Maven Central. Without a TornadoVM
runtime, or for small inputs, the classes run the lean Java path, which is also exact.

**Re-measured with TornadoVM 8.0.0** (`./bench.sh`, CUDA backend, same machine; 23/23 GPU parity tests pass on
CUDA and OpenCL): ASA 0.8 ms (1SMT) to 373 ms (3J3Q), 14.7-42.3x single-threaded BioJava and 2.1-3.1x lean Java;
all-pairs global 200 proteins 17.0 ms (586x BioJava, 24.5x lean Java), local 17.0 ms (615x, 17.2x); 500k pairs
up to 253 GCUPS. Pfam PF00104 global took 13.9 ms on 8.0.0 against 9 ms in the table above (7.0.1); not
investigated yet.

---

## Quick start

With [SDKMAN](https://sdkman.io), on Linux with an NVIDIA GPU (for AMD or Intel GPUs, use
`8.0.0-jdk21-opencl` and `TORNADO_BACKEND=opencl`):

```bash
curl -s "https://get.sdkman.io" | bash && source ~/.sdkman/bin/sdkman-init.sh   # if you have no SDKMAN yet
sdk install java 21.0.2-open
sdk install tornadovm 8.0.0-jdk21-cuda         # "n" when asked to make it the default is fine

git clone https://github.com/mikepapadim/biojava-tornadovm-extensions.git
cd biojava-tornadovm-extensions
./bench.sh            # builds, runs the GPU parity tests, then the ASA and alignment benchmarks
```

The scripts find TornadoVM 8.0.0 in SDKMAN by its version, whatever your default is, and install it
with `sdk install` if it is missing. Without SDKMAN they download the SDK from the GitHub release into
`.tornado/`. `tornado --devices` lists the devices; `-Dasa.fused.device=<driver>:<device>` and
`-Dsw.scores.device=<driver>:<device>` pick the device of the ASA and alignment kernels (default: TornadoVM's default device).

---

## For BioJava users

### 1. Add the library

It is not on Maven Central yet: install it locally (JDK 21 and Maven):

```bash
git clone https://github.com/mikepapadim/biojava-tornadovm-extensions.git
cd biojava-tornadovm-extensions && mvn -q install -DskipTests
```

```xml
<dependency>
  <groupId>org.biojava.tornado</groupId>
  <artifactId>biojava-tornadovm-extensions</artifactId>
  <version>0.1.0-SNAPSHOT</version>
</dependency>
```

It brings BioJava 7.3.0 (`biojava-structure`, `biojava-alignment`) with it. The TornadoVM API is a
`provided` dependency: at run time it comes from the TornadoVM SDK (step 3), so your application must
not bundle `tornado-api` itself.

### 2. Change one line

```java
// before: new AsaCalculator(structure, AsaCalculator.DEFAULT_PROBE_SIZE, 1000, nThreads, false)
TornadoAsaCalculator asa = new TornadoAsaCalculator(structure, AsaCalculator.DEFAULT_PROBE_SIZE,
        AsaCalculator.DEFAULT_N_SPHERE_POINTS, nThreads, false);
double[] perAtom = asa.calculateAsas();       // same as AsaCalculator.calculateAsas()
GroupAsa[] perResidue = asa.getGroupAsas();   // same as AsaCalculator.getGroupAsas()
```

The constructors mirror all of `AsaCalculator`'s: `Structure`, `Atom[]`, and `Point3d[]` with a fixed radius.

```java
// before: Alignments.getAllPairsScores(seqs, PairwiseSequenceScorerType.GLOBAL, gaps, matrix)
double[] scores = TornadoAlignments.getAllPairsScores(seqs, PairwiseSequenceScorerType.GLOBAL,
        new SimpleGapPenalty(), SubstitutionMatrixHelper.getBlosum62());

// before: Alignments.getAllPairsScorers(...) + Alignments.runPairwiseScorers(...)
List<PairwiseSequenceScorer<ProteinSequence, AminoAcidCompound>> scorers =
        TornadoAlignments.getAllPairsScorers(seqs, PairwiseSequenceScorerType.GLOBAL, gaps, matrix);
GuideTree<ProteinSequence, AminoAcidCompound> tree = new GuideTree<>(seqs, scorers);  // unchanged BioJava
```

Scores come in BioJava's order (pairs (i, j), i < j). The scorers carry the same score, max, min and
distance as BioJava's aligners, so anything built on them (guide trees, clustering) is unchanged. It works for
any compound set and substitution matrix (proteins, DNA with NUC.4.4, ambiguity codes).

### 3. Run with TornadoVM

Install TornadoVM 8.0.0 for JDK 21 with SDKMAN (`8.0.0-jdk21-cuda` for NVIDIA, `8.0.0-jdk21-opencl` for
AMD and Intel GPUs), or take the same SDK from the
[release](https://github.com/beehive-lab/TornadoVM/releases/tag/v8.0.0). Then start your application on
JDK 21 with the SDK's JVM options:

```bash
sdk install java 21.0.2-open
sdk install tornadovm 8.0.0-jdk21-cuda
export JAVA_HOME=$HOME/.sdkman/candidates/java/21.0.2-open
export TORNADOVM_HOME=$HOME/.sdkman/candidates/tornadovm/8.0.0-jdk21-cuda
# the SDK's options (module path, exports, native libraries), with its location filled in
sed "s|\${TORNADOVM_HOME}|$TORNADOVM_HOME|g" $TORNADOVM_HOME/tornado-argfile.template > tornado.args
$JAVA_HOME/bin/java @tornado.args -cp my-app.jar:<dependencies> com.example.Main
```

`scripts/run.sh` does exactly this for the classes of this project. The SDK puts its own JUnit (6.x)
on the module path, so tests run under TornadoVM use it rather than a second JUnit on the class path
(see `scripts/test-gpu.sh`).

Started with plain `java` (no argument file), the same program runs on the lean CPU path.

### When the GPU is used

| | ASA | All-pairs scores |
|---|---|---|
| GPU from | 1,000 atoms | 1e7 DP cells in total (about 250 pairs of 200 aa) |
| below that, or without TornadoVM | lean CPU path (`CpuAsa`), exact | lean CPU path (`CpuAlignmentScores`), exact |
| delegated to BioJava | - | `*_IDENTITIES`, `*_SIMILARITIES` (need a traceback), `KMERS`, `WU_MANBER` |

* `-Dbiojava.tornado=off` forces the CPU path, and `-Dbiojava.tornado=force` ignores the thresholds.
* GPUs without FP64 (e.g. Intel iGPUs) are detected: ASA then uses a single-precision kernel whose
  near-boundary cases are recomputed exactly on the CPU, so results stay identical.
* Tested on an RTX 4090 (CUDA and OpenCL backends) and an Intel UHD 770 (OpenCL); all tests are exact on all three.
* If a kernel fails (no device, out of memory, driver error), a warning is logged, the GPU path is
  disabled for the rest of the run, and the call is answered by the CPU path.
* The first GPU call in a JVM compiles the kernel (about 0.5 s). Later calls re-execute a cached execution
  plan, at about 1 ms (ASA) and 4 ms (alignment) of fixed overhead per call.
* Thread safety: calls are serialised on the device and are safe from several threads.

---

## For developers

```
src/main/java/org/biojava/tornado/
  TornadoSupport.java        GPU on/off, size thresholds, fallback
  asa/AsaKernels.java        fused kernel: cell-grid neighbour search, rank sort, occlusion test
  asa/TornadoAsaCalculator   host side, persistent capacity-padded execution plan
  align/SwKernels.java       strip-mined affine-gap NW/SW score kernel (one pair per thread)
  align/TornadoAlignments    encoding, cost-sorted batching, PrecomputedScorer
  asa/CpuAsa, align/CpuAlignmentScores   lean exact CPU paths (fallback, and the fair baseline)
  bench/                     AsaBench, SwBench
src/test/java/               parity tests: GPU result == BioJava result, exactly
```

Principles, which any new kernel should follow:
* **Exact parity.** The kernels reproduce BioJava's arithmetic: double precision in BioJava's
  operation order for ASA, and the same integer recurrences and boundary rules (`AlignerHelper`) for
  alignment. Tests compare with `assertArrayEquals(expected, actual, 0.0)`.
* **Drop-in API and fallback.** Same constructors and methods as the BioJava class; below a size threshold,
  without TornadoVM, and on any failure, run an exact CPU path.
* **Fair benchmarks.** Report the speedup against BioJava *and* against the same algorithm in lean Java
  on all cores.

```bash
./bench.sh                  # one command: TornadoVM SDK + data + build + GPU tests + benchmarks
./bench.sh asa | sw         # one benchmark
mvn test                    # plain JVM: exercises the CPU fallbacks
scripts/test-gpu.sh         # JUnit under TornadoVM: GPU results vs BioJava
scripts/run.sh <main> args  # run any class under TornadoVM
```

`scripts/env.sh` needs JDK 21 (`JAVA_HOME`, else SDKMAN's `21.0.2-open`) and uses `TORNADO_SDK` if set,
else SDKMAN's `tornadovm 8.0.0-jdk21-<backend>` (installed if missing), else downloads the release SDK
(`TORNADO_BACKEND=cuda|opencl`). Note that the `tornado` launcher prefers `TORNADOVM_HOME` over
`TORNADO_SDK`, so a stale `TORNADOVM_HOME` (for example from sdkman) silently runs another SDK; the scripts set both.

Hot spots investigated and not included are in [FINDINGS.md](FINDINGS.md) (CE structure alignment, and the
CUDA library bindings).

## Status and limitations

* A prototype proposed upstream in [biojava/biojava#1158](https://github.com/biojava/biojava/issues/1158).
  It is not on Maven Central.
* The alignment path covers score types `GLOBAL` and `LOCAL`, with affine, constant and linear gaps.
  `*_IDENTITIES` and `*_SIMILARITIES` (the MSA default) need a traceback and still run in BioJava.
* Local alignment on the GPU only matches parasail. Packed 8/16-bit arithmetic on the GPU would be the next step.
* CE structure alignment is not accelerated: its hot spot is `CECalculator.dpAlign`, see FINDINGS.md.
* Measured on one machine only (RTX 4090, Intel UHD 770). AMD and Apple GPUs have not been tested yet.
