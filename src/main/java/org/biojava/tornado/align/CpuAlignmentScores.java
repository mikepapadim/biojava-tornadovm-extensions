package org.biojava.tornado.align;

import java.util.stream.IntStream;

/**
 * Score-only affine-gap alignment on the CPU: the same recurrences as {@link SwKernels#alignmentScores} and as
 * BioJava's {@code AlignerHelper}, so the scores are exactly BioJava's, but in linear memory on primitive arrays
 * and without a traceback. Pairs run in parallel on all cores.
 */
public final class CpuAlignmentScores {

	private CpuAlignmentScores() {
	}

	/**
	 * @param seqs sequences as alphabet codes
	 * @param subs substitution scores, {@code subs[a*k+b]}
	 * @param k alphabet size
	 * @param gop gap open penalty (negative, as in BioJava)
	 * @param gep gap extension penalty (negative, as in BioJava)
	 * @param local Smith-Waterman if true, Needleman-Wunsch otherwise
	 * @return the score of every pair (i, j), i &lt; j, in BioJava's order
	 */
	public static int[] allPairs(int[][] seqs, int[] subs, int k, int gop, int gep, boolean local) {
		int n = seqs.length;
		int[] pairStart = new int[n];
		for (int i = 1; i < n; i++) {
			pairStart[i] = pairStart[i - 1] + (n - i);
		}
		int[] out = new int[n * (n - 1) / 2];
		IntStream.range(0, n).parallel().forEach(i -> {
			for (int j = i + 1; j < n; j++) {
				out[pairStart[i] + j - i - 1] = score(seqs[i], seqs[j], subs, k, gop, gep, local);
			}
		});
		return out;
	}

	/** The alignment score of query q against target t. */
	public static int score(int[] q, int[] t, int[] subs, int k, int gop, int gep, boolean local) {
		// as AlignerHelper: the lowest score that can still be extended without overflow
		int min = Integer.MIN_VALUE - gop - gep;
		int n = t.length;
		int[] maxRow = new int[n + 1], s0 = new int[n + 1], s1 = new int[n + 1];
		if (!local) {
			s1[0] = gop;
			int s2 = gop;
			for (int y = 1; y <= n; y++) {
				s2 += gep;
				s0[y] = min;
				s1[y] = min;
				maxRow[y] = s2;
			}
		}
		int best = 0;
		for (int x = 1; x <= q.length; x++) {
			int row = q[x - 1] * k;
			int diag = maxRow[0];
			int left0;
			int left2;
			if (local) {
				left0 = 0;
				left2 = 0;
				maxRow[0] = 0;
				s0[0] = 0;
				s1[0] = 0;
			} else {
				int v = s1[0] + gep;
				left0 = min;
				left2 = min;
				s0[0] = min;
				s1[0] = v;
				maxRow[0] = Math.max(v, min);
			}
			for (int y = 1; y <= n; y++) {
				int up0 = s0[y], up1 = s1[y], upM = maxRow[y];
				int a = diag + subs[row + t[y - 1]];
				int b = up1 >= up0 + gop ? up1 + gep : up0 + gop + gep;
				int c = left0 + gop >= left2 ? left0 + gop + gep : left2 + gep;
				if (local) {
					a = Math.max(a, 0);
					b = Math.max(b, 0);
					c = Math.max(c, 0);
					best = Math.max(best, a);
				}
				s0[y] = a;
				s1[y] = b;
				maxRow[y] = Math.max(Math.max(a, b), c);
				diag = upM;
				left0 = a;
				left2 = c;
			}
		}
		return local ? best : maxRow[n];
	}
}
