package org.biojava.tornado.asa;

import java.util.Arrays;
import java.util.stream.IntStream;

import javax.vecmath.Point3d;

/**
 * ASA on the CPU with the same algorithm as the fused GPU kernel (cell grid, per-atom neighbour list sorted by
 * distance, early-exit occlusion test) on primitive arrays. The per-atom results are bit-identical to
 * {@code AsaCalculator}'s: same neighbour criterion, same sphere points, same double-precision operation order.
 */
public final class CpuAsa {

	private CpuAsa() {
	}

	/**
	 * @param pts atom coordinates
	 * @param radii van der Waals radii
	 * @param probe probe radius
	 * @param nPoints number of sphere points
	 * @param parallel one task per atom on all cores if true, else single-threaded
	 * @return the ASA of every atom
	 */
	public static double[] calculate(Point3d[] pts, double[] radii, double probe, int nPoints, boolean parallel) {
		int n = pts.length;
		if (n == 0) {
			return new double[0];
		}
		double[] sphere = spherePoints(nPoints);
		double maxR = Arrays.stream(radii).max().orElse(0);
		double side = maxR + maxR + probe + probe;
		double lx = Double.MAX_VALUE, ly = lx, lz = lx, hx = -lx, hy = hx, hz = hx;
		double[] c = new double[3 * n];
		for (int i = 0; i < n; i++) {
			c[3 * i] = pts[i].x; c[3 * i + 1] = pts[i].y; c[3 * i + 2] = pts[i].z;
			lx = Math.min(lx, pts[i].x); ly = Math.min(ly, pts[i].y); lz = Math.min(lz, pts[i].z);
			hx = Math.max(hx, pts[i].x); hy = Math.max(hy, pts[i].y); hz = Math.max(hz, pts[i].z);
		}
		int nx = (int) ((hx - lx) / side) + 1, ny = (int) ((hy - ly) / side) + 1, nz = (int) ((hz - lz) / side) + 1;
		int[] cellOf = new int[n];
		int[] start = new int[nx * ny * nz + 1];
		for (int i = 0; i < n; i++) {
			cellOf[i] = ((int) ((c[3 * i] - lx) / side) * ny + (int) ((c[3 * i + 1] - ly) / side)) * nz + (int) ((c[3 * i + 2] - lz) / side);
			start[cellOf[i] + 1]++;
		}
		for (int k = 0; k + 1 < start.length; k++) start[k + 1] += start[k];
		int[] fill = Arrays.copyOf(start, start.length - 1);
		int[] sorted = new int[n];
		for (int i = 0; i < n; i++) sorted[fill[cellOf[i]]++] = i;
		final double fx = lx, fy = ly, fz = lz;
		final int fnx = nx, fny = ny, fnz = nz;
		double cons = 4.0 * Math.PI / nPoints;
		double[] asas = new double[n];
		IntStream atoms = IntStream.range(0, n);
		(parallel ? atoms.parallel() : atoms).forEach(i -> {
			double xi = c[3 * i], yi = c[3 * i + 1], zi = c[3 * i + 2];
			double ri = probe + radii[i];
			double[] nb = new double[4 * 64];
			double[] key = new double[64];
			int count = 0;
			int cx = (int) ((xi - fx) / side), cy = (int) ((yi - fy) / side), cz = (int) ((zi - fz) / side);
			for (int gx = cx - 1; gx <= cx + 1; gx++) for (int gy = cy - 1; gy <= cy + 1; gy++) for (int gz = cz - 1; gz <= cz + 1; gz++) {
				if (gx < 0 || gy < 0 || gz < 0 || gx >= fnx || gy >= fny || gz >= fnz) continue;
				int cell = (gx * fny + gy) * fnz + gz;
				for (int s = start[cell]; s < start[cell + 1]; s++) {
					int j = sorted[s];
					if (j == i) continue;
					double dx = c[3 * j] - xi, dy = c[3 * j + 1] - yi, dz = c[3 * j + 2] - zi;
					double d2 = dx * dx + dy * dy + dz * dz, d = Math.sqrt(d2);
					if (d < radii[i] + probe + probe + radii[j]) {
						if (count == key.length) { key = Arrays.copyOf(key, 2 * count); nb = Arrays.copyOf(nb, 8 * count); }
						double rj = radii[j] + probe;
						// insertion sort by distance
						int k = count - 1;
						while (k >= 0 && key[k] > d2) { key[k + 1] = key[k]; System.arraycopy(nb, 4 * k, nb, 4 * k + 4, 4); k--; }
						key[k + 1] = d2;
						nb[4 * k + 4] = dx; nb[4 * k + 5] = dy; nb[4 * k + 6] = dz;
						nb[4 * k + 7] = (d * d + ri * ri - rj * rj) / (2 * ri);
						count++;
					}
				}
			}
			int acc = 0;
			for (int p = 0; p < nPoints; p++) {
				double px = sphere[3 * p], py = sphere[3 * p + 1], pz = sphere[3 * p + 2];
				boolean ok = true;
				for (int o = 0; o < 4 * count; o += 4) {
					if (nb[o] * px + nb[o + 1] * py + nb[o + 2] * pz > nb[o + 3]) { ok = false; break; }
				}
				if (ok) acc++;
			}
			asas[i] = cons * acc * ri * ri;
		});
		return asas;
	}

	static double[] spherePoints(int nSpherePoints) {
		double[] points = new double[3 * nSpherePoints];
		double inc = Math.PI * (3.0 - Math.sqrt(5.0));
		double offset = 2.0 / nSpherePoints;
		for (int k = 0; k < nSpherePoints; k++) {
			double y = k * offset - 1.0 + (offset / 2.0);
			double r = Math.sqrt(1.0 - y * y);
			double phi = k * inc;
			points[3 * k] = Math.cos(phi) * r;
			points[3 * k + 1] = y;
			points[3 * k + 2] = Math.sin(phi) * r;
		}
		return points;
	}
}
