package org.biojava.tornado.bench;

import org.biojava.nbio.structure.Atom;
import org.biojava.nbio.structure.Structure;
import org.biojava.nbio.structure.StructureTools;
import org.biojava.nbio.structure.asa.AsaCalculator;
import org.biojava.tornado.asa.TornadoAsaCalculator;

/**
 * ASA parity check and timing: BioJava (1 thread, all cores) vs TornadoVM.
 * Usage: AsaBench PDBID [PDBID...]; property bench.reps (default 5).
 */
public class AsaBench {

	public static void main(String[] args) throws Exception {
		int reps = Integer.getInteger("bench.reps", 5);
		int cores = Runtime.getRuntime().availableProcessors();
		System.out.printf("%-6s %8s | %10s %10s | %10s %10s %10s | %9s %9s | %s%n", "pdb", "atoms", "cpu1 ms",
				"cpuN ms", "gpu ms", "(host)", "(tornado)", "vs cpu1", "vs cpuN", "parity");
		for (String id : args) {
			Structure s = Structures.load(id);
			Atom[] atoms = StructureTools.getAllNonHAtomArray(s, false, 0);

			double[] ref = null;
			double cpu1 = Double.MAX_VALUE, cpuN = Double.MAX_VALUE, gpu = Double.MAX_VALUE;
			for (int r = 0; r < reps; r++) {
				long t = System.nanoTime();
				ref = new AsaCalculator(atoms, AsaCalculator.DEFAULT_PROBE_SIZE, AsaCalculator.DEFAULT_N_SPHERE_POINTS, 1).calculateAsas();
				cpu1 = Math.min(cpu1, (System.nanoTime() - t) / 1e6);
			}
			for (int r = 0; r < reps; r++) {
				long t = System.nanoTime();
				new AsaCalculator(atoms, AsaCalculator.DEFAULT_PROBE_SIZE, AsaCalculator.DEFAULT_N_SPHERE_POINTS, cores).calculateAsas();
				cpuN = Math.min(cpuN, (System.nanoTime() - t) / 1e6);
			}
			double[] got = null;
			double nb = 0, tv = 0;
			for (int r = 0; r < reps + 1; r++) {
				TornadoAsaCalculator calc = new TornadoAsaCalculator(atoms, AsaCalculator.DEFAULT_PROBE_SIZE, AsaCalculator.DEFAULT_N_SPHERE_POINTS, cores);
				long t = System.nanoTime();
				got = calc.calculateAsasGpu();
				double ms = (System.nanoTime() - t) / 1e6;
				if (r == 0) {
					System.out.printf("  [%s cold GPU run incl. JIT: %.1f ms]%n", id, ms);
				} else if (ms < gpu) {
					gpu = ms;
					nb = calc.getLastHostNanos() / 1e6;
					tv = calc.getLastGpuNanos() / 1e6;
				}
			}
			System.out.printf("%-6s %8d | %10.1f %10.1f | %10.1f %10.1f %10.1f | %8.1fx %8.1fx | %s%n", id, atoms.length,
					cpu1, cpuN, gpu, nb, tv, cpu1 / gpu, cpuN / gpu, parity(ref, got));
		}
	}

	static String parity(double[] ref, double[] got) {
		double totRef = 0, totGot = 0, maxAbs = 0;
		int diff = 0;
		for (int i = 0; i < ref.length; i++) {
			totRef += ref[i];
			totGot += got[i];
			double d = Math.abs(ref[i] - got[i]);
			if (d > 1e-9) diff++;
			maxAbs = Math.max(maxAbs, d);
		}
		return String.format("total %.2f vs %.2f (rel %.1e), atoms differing %d, max |d| %.3f A^2", totRef, totGot,
				Math.abs(totRef - totGot) / totRef, diff, maxAbs);
	}
}
