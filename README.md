# biojava-tornadovm-extensions

GPU acceleration of two [BioJava](https://github.com/biojava/biojava) hot spots with
[TornadoVM](https://github.com/beehive-lab/TornadoVM), as **drop-in classes that return exactly the same
results** as the BioJava code they replace:

| BioJava class | Drop-in replacement | Result |
|---|---|---|
| `org.biojava.nbio.structure.asa.AsaCalculator` | `org.biojava.tornado.asa.TornadoAsaCalculator` | bit-identical per-atom and per-residue ASA |
| `org.biojava.nbio.alignment.Alignments.getAllPairsScores` / `getAllPairsScorers` (`GLOBAL`, `LOCAL`) | `org.biojava.tornado.align.TornadoAlignments` | identical scores, identical `GuideTree` |

On an RTX 4090 against a 32-thread i9-13900K (details and methodology in [FINDINGS.md](FINDINGS.md)):

| Workload | vs BioJava today | vs the same algorithm in lean Java on all 32 cores |
|---|---|---|
| ASA, 59k to 2.4M atoms | **11-17x** | 2.3-3.4x |
| All-pairs alignment scores, 20k to 2M pairs | **570-700x** (measured up to 40k pairs) | **15-56x** |
| Small inputs (< ~1k pairs) | the library runs its lean CPU code: 4-6x (ASA), ~25x (alignment) faster than BioJava | - |

BioJava is not modified. This library is built against BioJava 7.3.0 from Maven Central. Without a
TornadoVM runtime, or for small inputs, the classes run the same algorithms in lean Java on all cores. That
is also exact, and already several times faster than BioJava, so the same code runs everywhere.

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

Get a TornadoVM 7.0.1 SDK for JDK 21 from the
[releases](https://github.com/beehive-lab/TornadoVM/releases/tag/v7.0.1)
(`tornadovm-7.0.1-jdk21-cuda-linux-amd64` for NVIDIA, `-opencl-` for AMD/Intel, `-metal-mac-aarch64` for
Apple), or let `source scripts/env.sh` download it into `.tornado/`. Then start your application with
the SDK's argument file, on JDK 21:

```bash
export TORNADO_SDK=/path/to/tornadovm-7.0.1-cuda
java @$TORNADO_SDK/tornado-argfile -cp my-app.jar:<dependencies> com.example.Main
# or, equivalently: $TORNADO_SDK/bin/tornado -cp ... com.example.Main
```

Started with plain `java` (no argument file), the same program runs on the lean CPU path.

### When the GPU is used

| | ASA | All-pairs scores |
|---|---|---|
| GPU from | 1,000 atoms | 1e7 DP cells in total (about 250 pairs of 200 aa) |
| below that, or without TornadoVM | lean CPU path (`CpuAsa`), exact | lean CPU path (`CpuAlignmentScores`), exact |
| delegated to BioJava | - | `*_IDENTITIES`, `*_SIMILARITIES` (need a traceback), `KMERS`, `WU_MANBER` |

* `-Dbiojava.tornado=off` forces the CPU path, and `-Dbiojava.tornado=force` ignores the thresholds.
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

`scripts/env.sh` needs JDK 21 (`JAVA_HOME`) and uses `TORNADO_SDK` if set, else downloads the SDK
(`TORNADO_BACKEND=cuda|opencl`). Note that the `tornado` launcher prefers `TORNADOVM_HOME` over
`TORNADO_SDK`, so a stale `TORNADOVM_HOME` (for example from sdkman) silently runs another SDK; the scripts set both.

Hot spots investigated and not included are in [FINDINGS.md](FINDINGS.md) (CE structure alignment, and the
CUDA library bindings).
