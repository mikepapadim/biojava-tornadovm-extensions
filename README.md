# biojava-tornadovm-extensions

GPU acceleration of [BioJava](https://github.com/biojava/biojava) hot loops with
[TornadoVM](https://github.com/beehive-lab/TornadoVM), as drop-in classes that return
**exactly the same results** as the BioJava code they replace.

BioJava itself is not modified: this is a separate library, built against BioJava 7.3.0 from
Maven Central, with TornadoVM (JDK 21) as a `provided` dependency. Without a TornadoVM runtime,
or for small inputs, every class falls back to the stock BioJava implementation.

| BioJava | Drop-in | What runs on the GPU | Parity |
|---|---|---|---|
| `AsaCalculator` | `TornadoAsaCalculator` | neighbour search + Shrake-Rupley sphere-point occlusion, fused, one work-group per atom | bit-identical per-atom ASA (double precision, BioJava operation order) |
| `Alignments.getAllPairsScores` / `getAllPairsScorers` | `TornadoAlignments` | batched all-pairs Needleman-Wunsch / Smith-Waterman scores, affine gaps, one pair per thread, strip-mined | identical integer scores; identical `GuideTree` |

Measured speedups and methodology: see [FINDINGS.md](FINDINGS.md).

## Usage

```java
// ASA: same constructors and methods as AsaCalculator
double[] asas = new TornadoAsaCalculator(structure, AsaCalculator.DEFAULT_PROBE_SIZE,
        AsaCalculator.DEFAULT_N_SPHERE_POINTS, nThreads, false).calculateAsas();

// All-pairs scores (GLOBAL or LOCAL), same order as Alignments.getAllPairsScores
double[] scores = TornadoAlignments.getAllPairsScores(seqs, PairwiseSequenceScorerType.GLOBAL,
        new SimpleGapPenalty(), SubstitutionMatrixHelper.getBlosum62());

// ... or as scorers for a guide tree
GuideTree<ProteinSequence, AminoAcidCompound> tree = new GuideTree<>(seqs,
        TornadoAlignments.getAllPairsScorers(seqs, PairwiseSequenceScorerType.GLOBAL, gaps, blosum62));
```

## Quick start

Needs Linux x86-64, JDK 21, Maven and an NVIDIA GPU (or `TORNADO_BACKEND=opencl`).

```bash
./bench.sh          # downloads TornadoVM 7.0.1 into .tornado/, the PDB files, builds,
                    # runs the GPU parity tests and both benchmarks
./bench.sh asa      # or one benchmark: asa | sw
```

Any class of the project runs under TornadoVM with `scripts/run.sh <main-class> [args]`.
`scripts/env.sh` picks the JDK 21 and the SDK (set `TORNADO_SDK` to use an existing one).
Note that the `tornado` launcher prefers `TORNADOVM_HOME` over `TORNADO_SDK`: a stale
`TORNADOVM_HOME` (e.g. from sdkman) silently runs a different SDK.

## When the GPU is used

* Only under a TornadoVM runtime; `-Dbiojava.tornado=off` disables it, `-Dbiojava.tornado=force`
  ignores the size thresholds.
* Size thresholds (below them the CPU path is used): ASA 1000 atoms, all-pairs scores 5e7 DP cells in total.
* Not accelerated (always CPU): `*_IDENTITIES` /
  `*_SIMILARITIES` scorer types (they need a traceback); linear gap penalties.
* Any TornadoVM failure logs a warning, disables the GPU path for the rest of the run and falls
  back to BioJava.
* Each kernel is JIT-compiled once per JVM (~0.5 s) into a persistent execution plan that is
  re-executed for every later call; only the first call pays for it.

## Tests

```bash
mvn test              # plain JVM: exercises the CPU fallbacks
scripts/test-gpu.sh   # under TornadoVM: GPU results vs BioJava, exact equality
```
