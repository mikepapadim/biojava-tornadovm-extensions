package org.biojava.tornado.bench;

import java.util.stream.IntStream;

/**
 * A lean CPU baseline: the same score-only affine-gap DP as {@code SwKernels.alignmentScores} (and so BioJava's
 * scores), on primitive int arrays, one pair per task on all cores.
 */
final class CpuScores {

	private CpuScores() {
	}

	/** @param seqs sequences as alphabet codes; subs K*K; gop, gep negative */
	static int[] allPairs(int[][] seqs, int[] subs, int k, int gop, int gep, boolean local) {
		int n = seqs.length;
		int[] pairStart = new int[n];
		for (int i = 1; i < n; i++) pairStart[i] = pairStart[i - 1] + (n - i);
		int[] out = new int[n * (n - 1) / 2];
		IntStream.range(0, n).parallel().forEach(i -> {
			for (int j = i + 1; j < n; j++) out[pairStart[i] + j - i - 1] = score(seqs[i], seqs[j], subs, k, gop, gep, local);
		});
		return out;
	}

	static int score(int[] q, int[] t, int[] subs, int k, int gop, int gep, boolean local) {
		int min = Integer.MIN_VALUE - gop - gep;
		int n = t.length;
		int[] mArr = new int[n + 1], s0 = new int[n + 1], s1 = new int[n + 1];
		if (!local) {
			s1[0] = gop;
			int s2 = gop;
			for (int y = 1; y <= n; y++) { s2 += gep; s0[y] = min; s1[y] = min; mArr[y] = s2; }
		}
		int best = 0;
		for (int x = 1; x <= q.length; x++) {
			int row = q[x - 1] * k;
			int diag = mArr[0], left0, left2;
			if (local) { left0 = 0; left2 = 0; mArr[0] = 0; s0[0] = 0; s1[0] = 0; }
			else { int v = s1[0] + gep; left0 = min; left2 = min; s0[0] = min; s1[0] = v; mArr[0] = Math.max(v, min); }
			for (int y = 1; y <= n; y++) {
				int up0 = s0[y], up1 = s1[y], upM = mArr[y];
				int a = diag + subs[row + t[y - 1]];
				int b = up1 >= up0 + gop ? up1 + gep : up0 + gop + gep;
				int c = left0 + gop >= left2 ? left0 + gop + gep : left2 + gep;
				if (local) { a = Math.max(a, 0); b = Math.max(b, 0); c = Math.max(c, 0); best = Math.max(best, a); }
				s0[y] = a; s1[y] = b; mArr[y] = Math.max(Math.max(a, b), c);
				diag = upM; left0 = a; left2 = c;
			}
		}
		return local ? best : mArr[n];
	}
}
